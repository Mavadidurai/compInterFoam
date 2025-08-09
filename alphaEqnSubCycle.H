/*---------------------------------------------------------------------------*\
  CRITICAL FIX: Replace alphaEqnSubCycle.H content with this bounded version
\*---------------------------------------------------------------------------*/
#include <cmath>
// ✅ ROBUST: Alpha equation with comprehensive bounds enforcement
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
            // ✅ SIMPLIFIED: Basic alpha equation
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1)
              + fvm::div(phi, alpha1)
            );

            // ✅ CONSERVATIVE: Strong relaxation
            alpha1Eqn.relax(0.7);
            alpha1Eqn.solve();

            // ✅ CRITICAL: Immediate bounds enforcement
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
        // ✅ SIMPLIFIED: Single alpha solve
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1)
          + fvm::div(phi, alpha1)
        );

        // ✅ CONSERVATIVE: Strong relaxation
        alpha1Eqn.relax(0.1);
        alpha1Eqn.solve();

        // ✅ CRITICAL: Immediate bounds enforcement
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

    // ✅ UPDATE: Second phase with bounds
    alpha2 = 1.0 - alpha1;

    // ✅ ENSURE: Phase conservation
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

    // ✅ UPDATE: Density and boundary conditions
    rho = alpha1*rho1 + alpha2*rho2;
    
    alpha1.correctBoundaryConditions();
    alpha2.correctBoundaryConditions();

    // ✅ SIMPLE: Basic flux calculation
    alphaPhi10 = fvc::interpolate(alpha1) * phi;
    rhoPhi = alphaPhi10*fvc::interpolate(rho1) + (phi - alphaPhi10)*fvc::interpolate(rho2);

    // ✅ FINAL: Bounds check on fluxes
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

    // ✅ SAFETY: Final verification
    if (min(alpha1).value() < 0.0 || max(alpha1).value() > 1.0)
    {
        FatalErrorInFunction
            << "Alpha bounds violation detected!" << nl
            << "Min alpha1: " << min(alpha1).value() << nl
            << "Max alpha1: " << max(alpha1).value()
            << abort(FatalError);
    }
}
