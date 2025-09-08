#include <cmath>

{
    // Get time step info
    //const scalar deltaT = runTime.deltaT().value();
    
    if (verbose)
    {
        Info<< "Solving alpha equation with bounds enforcement" << endl;
    }

    // Runtime clamp controls
    const dictionary& solverControls =
        mesh.solutionDict().subOrEmptyDict("compInterFoamCoeffs");
    const dictionary& ctrl = runTime.controlDict();

    const Switch enableAlphaClamp = solverControls.lookupOrDefault<Switch>
    (
        "enableAlphaClamp",
        ctrl.lookupOrDefault<Switch>("enableAlphaClamp", true)
    );

    const scalar minAlpha = solverControls.lookupOrDefault<scalar>
    (
        "minAlpha",
        ctrl.lookupOrDefault<scalar>("minAlpha", 0.0)
    );

    const scalar maxAlpha = solverControls.lookupOrDefault<scalar>
    (
        "maxAlpha",
        ctrl.lookupOrDefault<scalar>("maxAlpha", 1.0)
    );

    label minClampCells = 0;
    label maxClampCells = 0;

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject("rhoPhiSum", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), Zero)
        );

        for (subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
             !(++alphaSubCycle).end();)
        {

            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1)
              + fvm::div(phi, alpha1)
            );


            alpha1Eqn.relax(0.7);
            alpha1Eqn.solve();


            if (enableAlphaClamp)
            {
                forAll(alpha1, cellI)
                {
                    if (!std::isfinite(alpha1[cellI]) || alpha1[cellI] < minAlpha)
                    {
                        alpha1[cellI] = minAlpha;
                        minClampCells++;
                    }
                    else if (alpha1[cellI] > maxAlpha)
                    {
                        alpha1[cellI] = maxAlpha;
                        maxClampCells++;
                    }
                }
            }
            else
            {
                forAll(alpha1, cellI)
                {
                    if (!std::isfinite(alpha1[cellI]))
                    {
                        alpha1[cellI] = minAlpha;
                    }
                }
            }

            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {

        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1)
          + fvm::div(phi, alpha1)
        );


        alpha1Eqn.relax(0.1);
        alpha1Eqn.solve();


        if (enableAlphaClamp)
        {
            forAll(alpha1, cellI)
            {
                if (!std::isfinite(alpha1[cellI]) || alpha1[cellI] < minAlpha)
                {
                    alpha1[cellI] = minAlpha;
                    minClampCells++;
                }
                else if (alpha1[cellI] > maxAlpha)
                {
                    alpha1[cellI] = maxAlpha;
                    maxClampCells++;
                }
            }
        }
        else
        {
            forAll(alpha1, cellI)
            {
                if (!std::isfinite(alpha1[cellI]))
                {
                    alpha1[cellI] = minAlpha;
                }
            }
        }
    }


    alpha2 = 1.0 - alpha1;


    if (enableAlphaClamp)
    {
        forAll(alpha2, cellI)
        {
            if (!std::isfinite(alpha2[cellI]) || alpha2[cellI] < minAlpha)
            {
                alpha2[cellI] = minAlpha;
                alpha1[cellI] = maxAlpha;
                maxClampCells++;
            }
            else if (alpha2[cellI] > maxAlpha)
            {
                alpha2[cellI] = maxAlpha;
                alpha1[cellI] = minAlpha;
                minClampCells++;
            }
        }
    }
    else
    {
        forAll(alpha2, cellI)
        {
            if (!std::isfinite(alpha2[cellI]))
            {
                alpha2[cellI] = maxAlpha;
                alpha1[cellI] = minAlpha;
            }
        }
    }

    reduce(minClampCells, sumOp<label>());
    reduce(maxClampCells, sumOp<label>());
    if (enableAlphaClamp && (minClampCells > 0 || maxClampCells > 0))
    {
        if (verbose)
        {
            Info<< "Alpha clamping stats: " << minClampCells
                << " cells to minAlpha (" << minAlpha << "), "
                << maxClampCells << " cells to maxAlpha (" << maxAlpha << ")" << endl;
        }
    }

    rho = alpha1*rho1 + alpha2*rho2;
    
    alpha1.correctBoundaryConditions();
    alpha2.correctBoundaryConditions();


    alphaPhi10 = fvc::interpolate(alpha1) * phi;
    rhoPhi = alphaPhi10*fvc::interpolate(rho1) + (phi - alphaPhi10)*fvc::interpolate(rho2);

    forAll(alphaPhi10, faceI)
    {
        if (!std::isfinite(alphaPhi10[faceI]))
        {
            alphaPhi10[faceI] = 0.0;
        }
    }
    if (verbose)
    {
        Info<< "Phase-1 volume fraction = " << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value() << endl;
    }
    if (min(alpha1).value() < 0.0 || max(alpha1).value() > 1.0)
    {
        WarningInFunction
            << "Alpha bounds violation detected after clamping:" << nl
            << "Min alpha1: " << min(alpha1).value() << nl
            << "Max alpha1: " << max(alpha1).value()
            << endl;
    }
}
