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
    const bool limitCompression = compressionLimiter < 1.0 - Foam::SMALL;
    
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

    scalar ocCoeff = 0;

    if
    (
        isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
     || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
    )
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

        if (verbose)
        {
            Info<< "alphaEqn: limiting interface compression with "
                << "compressionLimiter=" << compressionLimiter << endl;
        }
    }
    
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    if (scAlpha > 0)
    {
        phic += scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }

    const word sharpeningScheme
    (
        alphaControls.lookupOrDefault<word>("sharpeningScheme", "standard")
    );

    const scalar thinFilmCoeff
    (
        alphaControls.lookupOrDefault<scalar>("thinFilmCoeff", 0.0)
    );

    if (sharpeningScheme == "enhancedLimiter")
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

    surfaceScalarField phir(phic*nHatf);

    tmp<surfaceScalarField> phiAlpha(phiCN() + phir);
    const surfaceScalarField& phiAlphaRef = phiAlpha();

    #include "alphaSuSp.H"

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiAlpha,
                upwind<scalar>(mesh, phiAlpha)
            ).fvmDiv(phiAlphaRef, alpha1)
        );

        if (divU.valid())
        {
            alpha1Eqn -= Su() + fvm::Sp(Sp() + divU(), alpha1);
        }
        else
        {
            alpha1Eqn += fvm::Sp(Sp(), alpha1) + Su();
        }

        SolverPerformance<scalar> solverPerf = alpha1Eqn.solve(alpha1.name());

        alpha1 = Foam::max(Foam::min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        
        if (verbose)
        {
            Info<< "alpha1 solve: initial residual = "
                << solverPerf.initialResidual()
                << ", final residual = " << solverPerf.finalResidual()
                << ", iterations = " << solverPerf.nIterations() << endl;
        }

        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi10 = talphaPhi1UD();

        alpha2 = 1.0 - alpha1;
        alpha2.correctBoundaryConditions();
        
        mixture.correct();
    }

    for (label aCorr = 0; aCorr < nAlphaCorr; ++aCorr)
    {
        tmp<surfaceScalarField> talphaPhi1Un
        (
            fvc::flux
            (
                phiAlphaRef,
                (cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime())(),
                compressionFluxSchemeName
            )
        );

        if (MULESCorr)
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
                    (-Sp()*alpha1)(),
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

    #include "rhofs.H"

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
     || word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::localEulerDdtScheme<vector>::typeName
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

    rho = alpha1*rho1 + alpha2*rho2;

    if (verbose)
    {
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }
}
