/*---------------------------------------------------------------------------*\\
  Simplified alpha equation following the upstream compressibleInterFoam
  structure while keeping LIFT-specific hooks.
\\*---------------------------------------------------------------------------*/

extern Foam::Switch verbose;

{
    const dictionary& alphaControls = mesh.solverDict(alpha1.name());
    
    scalar compressionLimiter
    (
        alphaControls.lookupOrDefault<scalar>("compressionLimiter", 1.0)
    );

    const scalar minLimiter = Foam::SMALL;
    compressionLimiter = Foam::max(compressionLimiter, minLimiter);
    bool limitCompression = compressionLimiter < 1.0 - Foam::SMALL;
    
    #include "alphaScheme.H"

    tmp<fv::ddtScheme<scalar>> tDdtAlpha
    (
        fv::ddtScheme<scalar>::New
        (
            mesh,
            mesh.ddtScheme("ddt(alpha)")
        )
    );
    const fv::ddtScheme<scalar>& ddtAlpha = tDdtAlpha();
    const bool femtosecondDeltaT = runTime.deltaT().value() <= 1e-12;

    scalar ocCoeff = 0;

    if (isType<fv::EulerDdtScheme<scalar>>(ddtAlpha))
    {
        ocCoeff = 0;
    }
    else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
    {
        if (nAlphaSubCycles > 1)
        {
            FatalErrorInFunction
                << "Sub-cycling is not supported with the CrankNicolson ddt"
                << " scheme"
                << exit(FatalError);
        }

        if
        (
            alphaRestart
         || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
        )
        {
            ocCoeff =
                refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
               .ocCoeff();
        }
    }
    else
    {
        FatalErrorInFunction
            << "Only Euler and CrankNicolson ddt schemes are supported"
            << exit(FatalError);
    }

    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    tmp<surfaceScalarField> phiCN(phi);

    if (ocCoeff > 0 && phi.timeIndex() > 0)
    {
        phiCN = surfaceScalarField::New
        (
            "phiCN",
            cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime()
        );
    }

    surfaceScalarField phic
    (
        IOobject
        (
            "phic",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.cAlpha()*mag(phi/mesh.magSf())
    );
    
    tmp<surfaceScalarField> tMaxCompression;
    const surfaceScalarField* maxCompressionPtr = nullptr;

    if (limitCompression)
    {
        tMaxCompression = compressionLimiter*mag(phi)/mesh.magSf();
        maxCompressionPtr = &tMaxCompression();

        phic = Foam::min(phic, *maxCompressionPtr);

        if (verbose && Foam::Pstream::master())
        {
            Info<< "alphaEqn: limiting interface compression with "
                << "compressionLimiter=" << compressionLimiter << endl;
        }
    }
    // Per-call copies of the alpha controls. icAlpha, scAlpha and MULESCorr are
    // created once in createFields.H and live for the whole run; overwriting
    // them here made a single femtosecond step switch them off for every later
    // step (and a restart switch them on again).
    scalar icAlphaStep = icAlpha;
    scalar scAlphaStep = scAlpha;

    // The femtosecond branch used to also set compressionLimiter = 1 and
    // limitCompression = false. That came after the limiter had already been
    // applied above, so it never lifted the limit; it only skipped the second
    // (identical) application below. It is removed so the code does what the
    // documented and benchmarked settings say: effective compression
    // coefficient min(cAlpha, compressionLimiter) on every step.
    if (femtosecondDeltaT)
    {
        icAlphaStep = 0.0;
        scAlphaStep = 0.0;

        if (verbose && Foam::Pstream::master())
        {
            Info<< "alphaEqn: femtosecond timestep detected (dt = "
                << runTime.deltaT().value()
                << "); disabling secondary compression terms" << endl;
        }
    }
    if (icAlphaStep > 0)
    {
        phic *= (1.0 - icAlphaStep);
        phic += (mixture.cAlpha()*icAlphaStep)*fvc::interpolate(mag(U));
    }

    if (scAlphaStep > 0)
    {
        phic += scAlphaStep*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }

    const word sharpeningScheme
    (
        alphaControls.lookupOrDefault<word>("sharpeningScheme", "standard")
    );

    const scalar thinFilmCoeff
    (
        alphaControls.lookupOrDefault<scalar>("thinFilmCoeff", 0.0)
    );


    const bool useEnhancedLimiter =
        sharpeningScheme == "enhancedLimiter" && !femtosecondDeltaT;

    if (useEnhancedLimiter)
    {
        surfaceScalarField limiter
        (
            IOobject
            (
                "thinFilmLimiter",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            1.0 - thinFilmCoeff*fvc::interpolate(min(alpha1, alpha2))
        );

        limiter = Foam::max(Foam::min(limiter, scalar(1)), scalar(0));
        phic *= limiter;
    }
    else if (sharpeningScheme == "enhancedLimiter")
    {
        if (verbose && Foam::Pstream::master())
        {
            Info<< "alphaEqn: enhanced thin-film limiter disabled for femtosecond"
                " timestep to avoid over-sharpening" << endl;
        }
    }
    else if (sharpeningScheme != "standard")
    {
        WarningInFunction
            << "Unknown sharpeningScheme '" << sharpeningScheme
            << "' - using standard interface compression" << endl;
    }

    if (limitCompression)
    {
        phic = Foam::min(phic, *maxCompressionPtr);
    }

    // Option: restrict interface compression to the interface band.
    // The compression flux phic*nHatf*alpha*(1-alpha) is anti-diffusive
    // wherever alpha*(1-alpha) > 0, so inside the bulk film it sharpens
    // round-off-level alpha noise into growing dips (alpha 0.97-0.99)
    // that seed cell-to-cell velocity wiggles.  With the switch on,
    // compression acts only on faces touching a cell adjacent to an
    // alpha = 0.5 crossing (same band as the surface-tension band in pEqn.H).
    if (alphaControls.lookupOrDefault<Switch>("compressionInterfaceBand", true)   /* default on */)
    {
        volScalarField compressionBand
        (
            IOobject
            (
                "compressionBand",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        const scalarField& a1 = alpha1.primitiveField();
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        scalarField& band = compressionBand.primitiveFieldRef();

        forAll(nei, facei)
        {
            if ((a1[own[facei]] >= 0.5) != (a1[nei[facei]] >= 0.5))   // half-open: alpha = 0.5 counts once
            {
                band[own[facei]] = 1;
                band[nei[facei]] = 1;
            }
        }

        forAll(alpha1.boundaryField(), patchi)
        {
            const fvPatchScalarField& ap = alpha1.boundaryField()[patchi];
            if (ap.coupled())
            {
                const scalarField aNbr(ap.patchNeighbourField());
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                forAll(fc, i)
                {
                    if ((a1[fc[i]] >= 0.5) != (aNbr[i] >= 0.5))
                    {
                        band[fc[i]] = 1;
                    }
                }
            }
        }

        compressionBand.correctBoundaryConditions();
        phic *= localMax<scalar>(mesh).interpolate(compressionBand);
    }

    // use the mixture's interface normal flux, recomputed by
    // interfaceProperties::correct() at every alpha update. The solver used a
    // field built ONCE in createFields.H from the initial alpha, so the
    // compression acted only where the interface was at start-up/restart:
    // a moving film was never re-sharpened (cAlpha had no effect at all in a
    // translating-film test) and smeared into an alpha ~ 0.5 layer.
    surfaceScalarField phir(phic*mixture.nHatf());

    // the compression flux phir was added to the transport flux
    // (phiAlpha = phiCN + phir) for the upwind predictor, the high-order flux
    // and MULES::explicitSolve.  phir is not divergence-free, so alpha was
    // advected with a compressive velocity that does not vanish as alpha -> 1:
    // alpha overshoots 1 (to 1.5-2.7 in the oscillating-drop benchmark) and the
    // hard clips then delete metal (-9 % mass in 12 ns).  The stock form
    // transports alpha with phiCN only and adds the bounded alpha*(1 - alpha)
    // compression term to the high-order flux. Default on since it passed the
    // drop benchmarks (docs/VERIFICATION.md).
    const bool stockCompressionFlux =
        alphaControls.lookupOrDefault<Switch>("stockCompressionFlux", true);   // default on (bounded alpha transport)
    tmp<surfaceScalarField> phiAlpha
    (
        stockCompressionFlux ? tmp<surfaceScalarField>(phiCN()) : tmp<surfaceScalarField>(phiCN() + phir)
    );
    const surfaceScalarField& phiAlphaRef = phiAlpha();

    #include "alphaSuSp.H"

    // The MULES correction is skipped on femtosecond steps. The decision is
    // made per call and does not modify MULESCorr (see icAlphaStep above).
    const bool useMULESCorr = MULESCorr && !femtosecondDeltaT;
    const bool usedMULESCorrThisStep = useMULESCorr;

    if (!useMULESCorr)
    {
        // a cached correction from an earlier MULESCorr step must not be
        // re-applied after explicit steps in between
        talphaPhi1Corr0.clear();
    }

    if (useMULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiAlpha,
                upwind<scalar>(mesh, phiAlpha)
            ).fvmDiv(phiAlphaRef, alpha1)
        );

        if (divU.valid())
        {
            alpha1Eqn -= fvm::Sp(Sp() + divU(), alpha1);
            alpha1Eqn -= Su();
        }
        else
        {
            alpha1Eqn -= fvm::Sp(Sp(), alpha1);
            alpha1Eqn -= Su();
        }

        // NOTE: We intentionally avoid adding the historical "stabilizing" diagonal
        // term here.  On femtosecond timesteps (dt ~ 1e-14 s) the resulting
        // SMALL/dt contribution (O(1e-1)) overwhelms the physical diagonal (O(1e-3)),
        // forcing alpha1 to remain equal to its previous value and freezing the
        // phase-field evolution.
        
        SolverPerformance<scalar> solverPerf = alpha1Eqn.solve(alpha1.name());

        if (solverPerf.singular())
        {
            alpha1 = alpha1.oldTime();
            alpha1.correctBoundaryConditions();

            if (verbose && Foam::Pstream::master())
            {
                Info<< "alpha1 solve detected singular matrix; reverting to previous"
                    << " alpha field" << endl;
            }
        }

        alpha1 = Foam::max(Foam::min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        
        if (verbose && Foam::Pstream::master())
        {
            Info<< "alpha1 solve: initial residual = "
                << solverPerf.initialResidual()
                << ", final residual = " << solverPerf.finalResidual()
                << ", iterations = " << solverPerf.nIterations() << endl;
        }

        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi10 = talphaPhi1UD();

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            if (verbose && Foam::Pstream::master())
            {
                Info<< "Applying the previous iteration compression flux"
                    << endl;
            }

            if (divU.valid())
            {
                MULES::correct
                (
                    geometricOneField(),
                    alpha1,
                    talphaPhi1UD(),
                    talphaPhi1Corr0.ref(),
                    Sp(),
                    // The sources are already in the predictor; Su = -Sp*alpha1
                    // keeps them out of the correction (stock form). Passing Su()
                    // applied them again in every correction pass.
                    (-Sp()*alpha1())(),
                    oneField(),
                    zeroField()
                );
            }
            else
            {
                MULES::correct
                (
                    geometricOneField(),
                    alpha1,
                    talphaPhi1UD(),
                    talphaPhi1Corr0.ref(),
                    oneField(),
                    zeroField()
                );
            }

            alphaPhi10 += talphaPhi1Corr0();
        }

        // Cache this step's upwind flux; turned into the actual correction
        // (for re-use as "previous correction") once alphaPhi10 is final,
        // below.
        talphaPhi1Corr0 = talphaPhi1UD;

        alpha2 = 1.0 - alpha1;
        alpha2.correctBoundaryConditions();

        mixture.correct();
    }

    label effectiveAlphaCorr = nAlphaCorr;
    if (femtosecondDeltaT)
    {
        effectiveAlphaCorr = Foam::min(nAlphaCorr, label(1));
        if (MULESCorr && verbose && Foam::Pstream::master())
        {
            Info<< "alphaEqn: MULES correction skipped for this femtosecond"
                " timestep to reduce artificial compression" << endl;
        }
    }

    for (label aCorr = 0; aCorr < effectiveAlphaCorr; ++aCorr)
    {
        tmp<surfaceScalarField> talphaPhi1Un
        (
            stockCompressionFlux
          ? (
                fvc::flux
                (
                    phiCN(),
                    (cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime())(),
                    divAlphaName
                )
              + fvc::flux
                (
                    -fvc::flux(-phir, alpha2, "div(phirb,alpha)"),
                    alpha1,
                    "div(phirb,alpha)"
                )
            )
          : fvc::flux
            (
                phiAlphaRef,
                (cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime())(),
                compressionFluxSchemeName
            )
        );

        if (useMULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);
            volScalarField alpha10("alpha10", alpha1);

            if (divU.valid())
            {
                MULES::correct
                (
                    geometricOneField(),
                    alpha1,
                    talphaPhi1Un(),
                    talphaPhi1Corr.ref(),
                    Sp(),
                    (-Sp()*alpha1())(),   // sources already in the predictor
                    oneField(),
                    zeroField()
                );
            }
            else
            {
                MULES::correct
                (
                    geometricOneField(),
                    alpha1,
                    talphaPhi1Un(),
                    talphaPhi1Corr.ref(),
                    oneField(),
                    zeroField()
                );
            }

            if (aCorr == 0)
            {
                alphaPhi10 += talphaPhi1Corr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi10 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi10 = talphaPhi1Un();

            if (divU.valid())
            {
                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phiAlphaRef,
                    alphaPhi10,
                    Sp(),
                    (Su() + divU()*min(alpha1(), scalar(1)))(),
                    oneField(),
                    zeroField()
                );
            }
            else
            {
                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phiAlphaRef,
                    alphaPhi10,
                    Sp(),           // Include phase change diagonal
		    Su(),           // Include phase change source
                    oneField(),
                    zeroField()
                );
            }
        }

        alpha1 = Foam::max(Foam::min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();

        alpha2 = 1.0 - alpha1;
        alpha2.correctBoundaryConditions();
        
        mixture.correct();
    }

    if (alphaApplyPrevCorr && usedMULESCorrThisStep && talphaPhi1Corr0.valid())
    {
        // Turn the cached upwind flux into the actual correction applied
        // this step, so the next call seeds MULES with the real
        // compression correction rather than the raw upwind baseline.
        talphaPhi1Corr0 = alphaPhi10 - talphaPhi1Corr0;
    }

    #include "rhofs.H"

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN()*rho2f;
    }
    else
    {
        if (ocCoeff > 0)
        {
            alphaPhi10 =
                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
        }

        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
    }

    // rho is NOT updated here. Inside an alpha sub-cycle the
    // time index has advanced, so writing rho stored the intermediate
    // sub-cycle density as rho.oldTime(); ddt(rho,U) then no longer matched
    // div(rhoPhiSum) and a uniformly translating film carried a spurious
    // -rho_metal*U^2 pressure (Galilean test, cases/column1D_galilean_test). rho is set once
    // after the sub-cycles (compressibleAlphaEqnSubCycle.H), as in stock
    // compressibleInterFoam.

    if (verbose && Foam::Pstream::master())
    {
        Info<< "Phase-1 volume fraction (" << alpha1.name() << ") = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
        Info<< "Phase-2 volume fraction (" << alpha2.name() << ") = "
            << alpha2.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
            << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
            << endl;
    }
}
