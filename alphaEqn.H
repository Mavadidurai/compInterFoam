/*---------------------------------------------------------------------------*\
alphaEqn.H based on official OpenFOAM v2406 pattern
\*---------------------------------------------------------------------------*/
#include <cmath>
extern Foam::Switch verbose;
{
    word scheme("div(phi,alpha)");
    word schemer("div(phir,alpha)");

    tmp<fv::ddtScheme<scalar>> ddtAlpha
    (
        fv::ddtScheme<scalar>::New
        (
            mesh,
            mesh.ddtScheme("ddt(alpha)")
        )
    );

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    if
    (
        isType<fv::EulerDdtScheme<scalar>>(ddtAlpha())
     || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha())
    )
    {
        ocCoeff = 0;
    }
    else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha()))
    {
        if (nAlphaSubCycles > 1)
        {
            FatalErrorInFunction
                << "Sub-cycling is not supported "
                   "with the CrankNicolson ddt scheme"
                << exit(FatalError);
        }

        if
        (
            alphaRestart
         || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
        )
        {
            ocCoeff =
                refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha())
               .ocCoeff();
        }
    }
    else
    {
        FatalErrorInFunction
            << "Only Euler and CrankNicolson ddt schemes are supported"
            << exit(FatalError);
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);
        // Runtime selection of interface sharpening scheme
    const dictionary& alphaControls = mesh.solverDict(alpha1.name());
    const word sharpeningScheme
    (
        alphaControls.lookupOrDefault<word>("sharpeningScheme", "standard")
    );
    const scalar thinFilmCoeff
    (
        alphaControls.lookupOrDefault<scalar>("thinFilmCoeff", 0.0)
    );
    const scalar alphaSourceFloor
    (
        alphaControls.lookupOrDefault<scalar>("alphaSourceFloor", 1e-3) // recommended ~1e-3
    );
    const scalar sourceMaxMag
    (
        alphaControls.lookupOrDefault<scalar>("sourceMaxMag", GREAT)
    );
    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Add the optional shear compression contribution
    if (scAlpha > 0)
    {
        phic +=
            scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }

    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    forAll(phicBf, patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];
        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    // Use surfaceScalarField instead of tmp<>
    surfaceScalarField phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0 && phi.timeIndex() > 0)
    {
        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
    }

    // Construct interface compression flux with optional limiter
    surfaceScalarField phir(phic*nHatf);

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

        limiter = max(min(limiter, scalar(1)), scalar(0));
        phir = phic*limiter*nHatf;
    }
    else if (sharpeningScheme == "isoAdvector")
    {
        Info<< "isoAdvector selected but not available - using standard scheme" << endl;
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        volScalarField alpha10(alpha1);

        // Remove complex source terms for stability
        volScalarField::Internal SpInternal
        (
            IOobject
            (
                "SpInternal",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dgdt.dimensions(), Zero)
        );

        volScalarField::Internal SuInternal
        (
            IOobject
            (
                "SuInternal",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dgdt.dimensions(), Zero)
        );

               // Basic source term calculation with clamping
        label sourceClampCells = 0;
        bool hitAlphaFloor = false;
        scalarField& SpInternalField = SpInternal.primitiveFieldRef();
        scalarField& SuInternalField = SuInternal.primitiveFieldRef();
        const scalarField& alpha1Field = alpha1.internalField();
        const scalarField& dgdtField = dgdt.internalField();

        forAll(dgdtField, celli)
        {
            const scalar dgt = dgdtField[celli];

            // Separate evaporation (positive) and condensation (negative) contributions
            scalar evapRate = max(dgt, scalar(0));   // posPart(dgdt)
            scalar condRate = max(-dgt, scalar(0));  // posPart(-dgdt) ~ negPart(dgdt)

            bool cellClamped = false;
            if (evapRate > sourceMaxMag)
            {
                evapRate = sourceMaxMag;
                cellClamped = true;
            }
            if (condRate > sourceMaxMag)
            {
                condRate = sourceMaxMag;
                cellClamped = true;
            }

            scalar evapDenom = max(alpha1Field[celli], alphaSourceFloor);
            if (evapRate > 0 && alpha1Field[celli] < alphaSourceFloor)
            {
                hitAlphaFloor = true;
            }
            if (evapDenom > sourceMaxMag)
            {
                evapDenom = sourceMaxMag;
                cellClamped = true;
            }

            const scalar alpha2Cell = 1.0 - alpha1Field[celli];
            scalar condDenom = max(alpha2Cell, alphaSourceFloor);
            if (condRate > 0 && alpha2Cell < alphaSourceFloor)
            {
                hitAlphaFloor = true;
            }
            if (condDenom > sourceMaxMag)
            {
                condDenom = sourceMaxMag;
                cellClamped = true;
            }

            const scalar evapCoeff = evapRate/evapDenom;
            const scalar condCoeff = condRate/condDenom;

            // Implicit Sp terms retain the true evaporation/condensation scaling
            // and reduce to the previous dgt/denom formulation
            SpInternalField[celli] = -(evapCoeff + condCoeff);
            SuInternalField[celli] = condCoeff;

            if (cellClamped)
            {
                sourceClampCells++;
            }
        }
        reduce(sourceClampCells, sumOp<label>());
        if (sourceClampCells > 0)
        {
            WarningInFunction
                << "Clamped " << sourceClampCells
                << " source term values to max magnitude "
                << sourceMaxMag << endl;
        }
        if (hitAlphaFloor)
        {
            WarningInFunction
                << "alphaSourceFloor triggered in source term division;"
                << " solution may be stiff" << endl;
        }

        label nonFiniteSources = 0;
        forAll(SpInternal, celli)
        {
            if (!std::isfinite(SpInternal[celli]) || !std::isfinite(SuInternal[celli]))
            {
                SpInternal[celli] = 0.0;
                SuInternal[celli] = 0.0;
                nonFiniteSources++;
            }
        }
        reduce(nonFiniteSources, sumOp<label>());
        if (nonFiniteSources > 0)
        {
            WarningInFunction
                << nonFiniteSources
                << " non-finite source term values set to zero" << endl;
        }
        volScalarField spField
        (
            IOobject(
                "spField",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            SpInternal,
            alpha1.boundaryField()
        );

        volScalarField suField
        (
            IOobject(
                "suField",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            SuInternal,
            alpha1.boundaryField()
        );
        // Create proper volScalarField for divU
        volScalarField divUField
        (
            IOobject
            (
                "divUField",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh.moving()
          ? fvc::div(phiCN + mesh.phi())
          : fvc::div(phiCN)
        );

        // Create surface fields for MULES
        surfaceScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        );

        surfaceScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        );

        Su = fvc::interpolate
        (
            divUField*min(alpha1, dimensionedScalar("one", dimless, 1.0))
        );
        
        Sp += fvc::interpolate(spField);
        Su += fvc::interpolate(suField);

        // Basic alpha equation
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
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
        - fvm::div(phir, alpha1)
        + fvm::Sp(SpInternal, alpha1) + SuInternal
        );
               alpha1Eqn += fvm::Sp
        (
            dimensionedScalar("diagStab", alpha1.dimensions()/dimTime, VSMALL),
            alpha1
        );

        // Add small stabilization term to diagonal
        alpha1Eqn += fvm::Sp(dimensionedScalar("stab", dimless/dimTime, 1e-9), alpha1);
        
        // Solve with residual checking
        SolverPerformance<scalar> solverPerf = alpha1Eqn.solve(alpha1.name());
        
        if (solverPerf.initialResidual() > 1e-4 && !solverPerf.converged())
        {
            WarningInFunction
                << "Alpha equation may not have converged"
                << "  Initial residual: " << solverPerf.initialResidual()
                << "  Final residual: " << solverPerf.finalResidual()
                << "  Iterations: " << solverPerf.nIterations()
                << endl;
        }
        
        // Check for non-finite values before bounding
        {
            const scalarField& alpha1Field = alpha1.internalField();
            bool hasNonFinite = false;
            forAll(alpha1Field, i)
            {
                if (!std::isfinite(alpha1Field[i]))
                {
                    hasNonFinite = true;
                    break;
                }
            }
            if (hasNonFinite)
            {
                FatalErrorInFunction
                    << "Non-finite values detected in alpha1 field after solve"
                    << exit(FatalError);
            }
        }

        // Bound and update fields with validation
        alpha1 = min(max(alpha1, scalar(0)), scalar(1));
        alpha2 = 1.0 - alpha1;
        
        // Validate density fields
        if (rho1.size() != mesh.nCells() || rho2.size() != mesh.nCells())
        {
            FatalErrorInFunction
                << "Density field size mismatch: rho1:" << rho1.size()
                << " rho2:" << rho2.size()
                << " mesh:" << mesh.nCells()
                << exit(FatalError);
        }

        // Check for non-finite density values
        const scalarField& rho1Internal = rho1.internalField();
        const scalarField& rho2Internal = rho2.internalField();
        forAll(rho1Internal, cellI)
        {
            if (!std::isfinite(rho1Internal[cellI]) || !std::isfinite(rho2Internal[cellI]))
            {
                FatalErrorInFunction
                    << "Non-finite density values detected at cell " << cellI
                    << " rho1:" << rho1Internal[cellI]
                    << " rho2:" << rho2Internal[cellI]
                    << exit(FatalError);
            }
        }
        
        rho = alpha1*rho1 + alpha2*rho2;
        alpha1.correctBoundaryConditions();
        alpha2.correctBoundaryConditions();

        // Get flux safely without tmp issues
        alphaPhi10 = alpha1Eqn.flux();

        if (MULESCorr)
        {
            // Optionally revert to previous iteration values
            if (alphaApplyPrevCorr)
            {
                alpha1 = alpha10;
            }

            surfaceScalarField alphaPhiCorr(alphaPhi10);

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                phiCN + phir,
                alphaPhiCorr,
                Sp,
                Su,
                oneField(),
                zeroField()
            );

            alphaPhi10 = alphaPhiCorr;
        }
        else
        {
            // Basic MULES explicit solve without correction
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phiCN + phir,
                alphaPhi10,
                Sp,
                Su,
                oneField(),
                zeroField()
            );
        }
        alpha2 = 1.0 - alpha1;

        if (verbose)
        {
            Info<< "Phase-1 volume fraction = "
                << alpha1.weightedAverage(mesh.Vsc()).value()
                << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
                << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
                << endl;
        }
    }

    tmp<surfaceScalarField> trho1 = fvc::interpolate(rho1);
    tmp<surfaceScalarField> trho2 = fvc::interpolate(rho2);
    rhoPhi = alphaPhi10*(trho1() - trho2()) + phiCN*trho2();

    if (ocCoeff > 0 && phi.timeIndex() > 0)
    {
        // Calculate the end-of-time-step alpha flux
        alphaPhi10 = (phi - (1.0 - cnCoeff)*phi.oldTime())/cnCoeff;
    }
}

