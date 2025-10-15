extern Foam::Switch verbose;
extern const bool master;

alpha1 = max(min(alpha1, scalar(1)), scalar(0));
alpha1.correctBoundaryConditions();
alpha2 = scalar(1) - alpha1;
alpha2.correctBoundaryConditions();

const label electronSubCyclesForAlpha =
    // Use the const helper so we base the guard on the pending global step
    ttm.plannedElectronSubCycles(runTime.deltaT());
const label nAlphaSubCyclesForStep =
    Foam::max(nAlphaSubCycles, electronSubCyclesForAlpha);

if (useAdvancedCapturing && pInterfaceCapturing.valid())
{
    tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

    if (nAlphaSubCyclesForStep > 1)
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
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCyclesForStep);
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

    if (nAlphaSubCyclesForStep <= 1)
    {
        tmp<surfaceScalarField> trho1 = fvc::interpolate(rho1);
        tmp<surfaceScalarField> trho2 = fvc::interpolate(rho2);
        rhoPhi = alphaPhi10*trho1() + (phi - alphaPhi10)*trho2();
    }

if (verbose && master)
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
else
{
    tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

    if (nAlphaSubCyclesForStep > 1)
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
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCyclesForStep);
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
