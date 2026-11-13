extern Foam::Switch verbose;

// Alpha sub-cycling is controlled by fvSolution (nAlphaSubCycles) alone;
// the electron micro-step count is a thermal-stiffness control and has no
// bearing on the interface CFL, so tying the two multiplied the cost of
// the alpha transport (and the per-sub-cycle recoil updates) needlessly.
const label nAlphaSubCyclesForStep = nAlphaSubCycles;

// One code path for both settings. The two branches used to be ~90-line
// copies that differed only in the interface-capturing update, the rhoPhi
// reconstruction without sub-cycles and the verbose report; those three
// steps are kept, and still apply only with advanced capturing.
const bool advancedCapturingActive =
    useAdvancedCapturing && pInterfaceCapturing.valid();

{
    tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

    // bound alpha1 to [0, 1] and refresh alpha2
    auto boundAlpha = [&]()
    {
        alpha1 = max(min(alpha1, scalar(1)), scalar(0));
        alpha1.correctBoundaryConditions();
        alpha2 = scalar(1) - alpha1;
        alpha2.correctBoundaryConditions();
    };

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
            if (advancedCapturingActive)
            {
                pInterfaceCapturing->correct();
            }

            // alpha1 is bounded once before the transport; alphaEqn.H bounds
            // it again after each solve, so no further clip is needed here
            boundAlpha();
            #include "alphaEqn.H"

            const dimensionedScalar subCycleWeight(runTime.deltaT()/totalDeltaT);
            talphaPhi1.ref() += subCycleWeight*alphaPhi10;
            rhoPhiSum += subCycleWeight*rhoPhi;
        }

        alphaPhi10 = talphaPhi1();
        rhoPhi = rhoPhiSum;
    }
    else
    {
        if (advancedCapturingActive)
        {
            pInterfaceCapturing->correct();
        }

        boundAlpha();
        #include "alphaEqn.H"
    }

    rho = alpha1*rho1 + alpha2*rho2;

    if (advancedCapturingActive && nAlphaSubCyclesForStep <= 1)
    {
        tmp<surfaceScalarField> trho1 = fvc::interpolate(rho1);
        tmp<surfaceScalarField> trho2 = fvc::interpolate(rho2);
        rhoPhi = alphaPhi10*trho1() + (phi - alphaPhi10)*trho2();
    }

    if (advancedCapturingActive && verbose && Foam::Pstream::master())
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
