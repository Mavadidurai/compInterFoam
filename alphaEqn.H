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
        const scalar diagStab
    (
        alphaControls.lookupOrDefault<scalar>("diagStab", SMALL)
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
        forAll(dgdt, celli)
        {
            scalar dgt = dgdt[celli];
            if (dgt > 0.0)
            {
                if (alpha1[celli] < alphaSourceFloor)
                {
                    hitAlphaFloor = true;
                }
                scalar denom = max(alpha1[celli], alphaSourceFloor);
                bool cellClamped = false;
                if (mag(dgt) > sourceMaxMag)
                {
                    dgt = sign(dgt)*sourceMaxMag;
                    cellClamped = true;
                }
                if (denom > sourceMaxMag)
                {
                    denom = sourceMaxMag;
                    cellClamped = true;
                }
                if (cellClamped)
                {
                    sourceClampCells++;
                }
                // Evaporation: remove mass from liquid
                SpInternal[celli] -= dgt/denom;
            }
            else if (dgt < 0.0)
            {
                scalar denom = max(1.0 - alpha1[celli], alphaSourceFloor);
                if ((1.0 - alpha1[celli]) < alphaSourceFloor)
                {
                    hitAlphaFloor = true;
                }
                bool cellClamped = false;
                if (mag(dgt) > sourceMaxMag)
                {
                    dgt = sign(dgt)*sourceMaxMag;
                    cellClamped = true;
                }
                if (denom > sourceMaxMag)
                {
                    denom = sourceMaxMag;
                    cellClamped = true;
                }
                if (cellClamped)
                {
                    sourceClampCells++;
                }
                // Condensation: add mass to liquid
                SpInternal[celli] += dgt/denom;
                SuInternal[celli] -= dgt/denom;
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

        Su = fvc::interpolate(divUField*min(alpha1, scalar(1)));
        
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
        // Ensure a strictly positive diagonal with proper dimensions
        alpha1Eqn += fvm::Sp(
            dimensionedScalar("stabCoeff", alpha1.dimensions()/dimTime, diagStab),
            alpha1
        );
        
        if (Foam::gMin(alpha1Eqn.diag()) <= 0)
        {
            Warning
                << "alpha1Eqn: detected non-positive diagonal entry."
                << " Adjusting to SMALL." << nl;
            scalarField& diag = alpha1Eqn.diag();
            forAll(diag, i)
            {
                if (diag[i] <= 0)
                {
                    diag[i] = SMALL;
                }
            }
        }

        alpha1Eqn.solve();

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

    rhoPhi = alphaPhi10*(fvc::interpolate(rho1) - fvc::interpolate(rho2)) + phiCN*fvc::interpolate(rho2);

    if (ocCoeff > 0 && phi.timeIndex() > 0)
    {
        // Calculate the end-of-time-step alpha flux
        alphaPhi10 = (phi - (1.0 - cnCoeff)*phi.oldTime())/cnCoeff;
    }
}
