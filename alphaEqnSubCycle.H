#include <cmath>

{
    // Get time step info
    //const scalar deltaT = runTime.deltaT().value();
    
    Info<< "Solving alpha equation with bounds enforcement" << endl;
    
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


            forAll(alpha1, cellI)
            {
                if (alpha1[cellI] < 0.0) 
                {
                    alpha1[cellI] = 0.0;
                }
                else if (alpha1[cellI] > 1.0) 
                {
                    alpha1[cellI] = 1.0;
                }
                else if (!std::isfinite(alpha1[cellI]))
                {
                    alpha1[cellI] = 0.0;
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


        forAll(alpha1, cellI)
        {
            if (alpha1[cellI] < 0.0) 
            {
                alpha1[cellI] = 0.0;
            }
            else if (alpha1[cellI] > 1.0) 
            {
                alpha1[cellI] = 1.0;
            }
            else if (!std::isfinite(alpha1[cellI]))
            {
                alpha1[cellI] = 0.0;
            }
        }
    }


    alpha2 = 1.0 - alpha1;


    forAll(alpha2, cellI)
    {
        if (alpha2[cellI] < 0.0) 
        {
            alpha2[cellI] = 0.0;
            alpha1[cellI] = 1.0;
        }
        else if (alpha2[cellI] > 1.0) 
        {
            alpha2[cellI] = 1.0;
            alpha1[cellI] = 0.0;
        }
        else if (!std::isfinite(alpha2[cellI]))
        {
            alpha2[cellI] = 1.0;
            alpha1[cellI] = 0.0;
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

    Info<< "Phase-1 volume fraction = " << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value() << endl;


    if (min(alpha1).value() < 0.0 || max(alpha1).value() > 1.0)
    {
        WarningInFunction
            << "Alpha bounds violation detected after clamping:" << nl
            << "Min alpha1: " << min(alpha1).value() << nl
            << "Max alpha1: " << max(alpha1).value()
            << endl;
    }
}
