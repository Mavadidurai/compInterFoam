/*---------------------------------------------------------------------------*\
alphaEqn.H based on official OpenFOAM v2406 pattern
\*---------------------------------------------------------------------------*/

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

    // Standard face-flux compression coefficient
    surfaceScalarField phic(antiDiffusion*mixture.cAlpha()*mag(phi/mesh.magSf()));

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

    surfaceScalarField phir(phic*nHatf);

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

        // Basic source term calculation
        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0)
            {
                SpInternal[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
                SuInternal[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
            }
            else if (dgdt[celli] < 0.0)
            {
                SpInternal[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
            }
        }

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

        // Basic alpha equation without complex MULES correction
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
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        // Get flux safely without tmp issues
        alphaPhi10 = alpha1Eqn.flux();

        // Basic MULES explicit solve without correction
        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1,
            phiCN,
            alphaPhi10,
            Sp,
            Su,
            oneField(),
            zeroField()
        );

        alpha2 = 1.0 - alpha1;

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }

    rhoPhi = alphaPhi10*(fvc::interpolate(rho1) - fvc::interpolate(rho2)) + phiCN*fvc::interpolate(rho2);

    if (ocCoeff > 0 && phi.timeIndex() > 0)
    {
        // Calculate the end-of-time-step alpha flux
        alphaPhi10 = (phi - (1.0 - cnCoeff)*phi.oldTime())/cnCoeff;
    }
}
