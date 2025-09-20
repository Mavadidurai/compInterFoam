alpha1 = max(min(alpha1, scalar(1)), scalar(0));
alpha1.correctBoundaryConditions();
alpha2 = scalar(1) - alpha1;
alpha2.correctBoundaryConditions();

if (pInterfaceCapturing.valid())
{
    tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

    if (nAlphaSubCycles > 1)
    {
        const dimensionedScalar totalDeltaT(runTime.deltaT());

        talphaPhi1.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "alphaPhi1",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(alphaPhi10.dimensions(), Zero)
            )
        );

        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), Zero)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            pInterfaceCapturing->correct();

            alpha1 = max(min(alpha1, scalar(1)), scalar(0));
            alpha1.correctBoundaryConditions();
            alpha2 = scalar(1) - alpha1;
            alpha2.correctBoundaryConditions();

            #include "alphaEqn.H"

            alpha1 = max(min(alpha1, scalar(1)), scalar(0));
            alpha1.correctBoundaryConditions();
            alpha2 = scalar(1) - alpha1;
            alpha2.correctBoundaryConditions();

            const dimensionedScalar subCycleWeight(runTime.deltaT()/totalDeltaT);
            talphaPhi1.ref() += subCycleWeight*alphaPhi10;
            rhoPhiSum += subCycleWeight*rhoPhi;
        }

        alphaPhi10 = talphaPhi1();
        rhoPhi = rhoPhiSum;
    }
    else
    {
        pInterfaceCapturing->correct();

        alpha1 = max(min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        alpha2 = scalar(1) - alpha1;
        alpha2.correctBoundaryConditions();

        #include "alphaEqn.H"

        alpha1 = max(min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        alpha2 = scalar(1) - alpha1;
        alpha2.correctBoundaryConditions();
    }

    rho = alpha1*rho1 + alpha2*rho2;

    if (nAlphaSubCycles <= 1)
    {
        tmp<surfaceScalarField> trho1 = fvc::interpolate(rho1);
        tmp<surfaceScalarField> trho2 = fvc::interpolate(rho2);
        rhoPhi = alphaPhi10*trho1() + (phi - alphaPhi10)*trho2();
    }

    if (verbose)
    {
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value() << endl;
    }
}
else
{
    tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

    if (nAlphaSubCycles > 1)
    {
        const dimensionedScalar totalDeltaT(runTime.deltaT());

        talphaPhi1.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "alphaPhi1",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(alphaPhi10.dimensions(), Zero)
            )
        );

        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), Zero)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            alpha1 = max(min(alpha1, scalar(1)), scalar(0));
            alpha1.correctBoundaryConditions();
            alpha2 = scalar(1) - alpha1;
            alpha2.correctBoundaryConditions();

            #include "alphaEqn.H"

            alpha1 = max(min(alpha1, scalar(1)), scalar(0));
            alpha1.correctBoundaryConditions();
            alpha2 = scalar(1) - alpha1;
            alpha2.correctBoundaryConditions();

            const dimensionedScalar subCycleWeight(runTime.deltaT()/totalDeltaT);
            talphaPhi1.ref() += subCycleWeight*alphaPhi10;
            rhoPhiSum += subCycleWeight*rhoPhi;
        }

        alphaPhi10 = talphaPhi1();
        rhoPhi = rhoPhiSum;
    }
    else
    {
        alpha1 = max(min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        alpha2 = scalar(1) - alpha1;
        alpha2.correctBoundaryConditions();

        #include "alphaEqn.H"

        alpha1 = max(min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        alpha2 = scalar(1) - alpha1;
        alpha2.correctBoundaryConditions();
    }

    rho = alpha1*rho1 + alpha2*rho2;
}
