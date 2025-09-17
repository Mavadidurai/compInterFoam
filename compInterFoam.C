/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) OpenCFD OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    compressibleInterFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Either mixture or two-phase transport modelling may be selected.  In the
    mixture approach a single laminar, RAS or LES model is selected to model the
    momentum stress.  In the Euler-Euler two-phase approach separate laminar,
    RAS or LES selected models are selected for each of the phases.

    Extended with femtosecond laser modeling, two-temperature physics, 
    and advanced interface capturing for Laser-Induced Forward Transfer (LIFT).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "autoPtr.H"


#include "fvcSmooth.H"
#include "twoPhaseMixtureThermo.H"
#include "extrapolatedCalculatedFvPatchFields.H"
// LIFT-specific model headers
#include "femtosecondLaserModel.H"
#include "twoTemperatureModel.H"
#include "advancedInterfaceCapturing.H"

Foam::Switch verbose(false);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two compressible, non-isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing."
        " Extended with LIFT (Laser-Induced Forward Transfer) models."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"
    #ifndef CREATE_FIELDS_DONE
    #include "createFields.H"
    #define CREATE_FIELDS_DONE
    #endif
    
    // Track recoil pressure update intervals when alpha subcycle is skipped
    label recoilCallCount = 0;

    // Reference to psi fields (needed for compressibility)
    #ifndef NDEBUG
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();
    // Silence unused variable warnings in debug mode
    (void)psi1;
    (void)psi2;
    #endif

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
    // Instantiate PIMPLE control for pressure-velocity coupling
    pimpleControl pimple(mesh);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            // Traditional inline LTS implementation
            volScalarField& rDeltaT = trDeltaT.ref();
            
            const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
            
            scalar maxCo = pimpleDict.getOrDefault<scalar>("maxCo", 0.9);
            scalar maxAlphaCo = pimpleDict.getOrDefault<scalar>("maxAlphaCo", 0.2);
            
            // Set reciprocal time-step from local Courant number
            rDeltaT.ref() = max
            (
                1/dimensionedScalar("maxDeltaT", dimTime, GREAT),
                fvc::surfaceSum(mag(rhoPhi))()()
               /((2*maxCo)*mesh.V()*rho())
            );
            
            // Limit based on alpha Courant number
            if (maxAlphaCo < maxCo)
            {
                volScalarField alpha1Bar(fvc::average(alpha1));
                rDeltaT.ref() = max
                (
                    rDeltaT(),
                    pos0(alpha1Bar() - 0.01)
                   *pos0(0.99 - alpha1Bar())
                   *fvc::surfaceSum(mag(phi))()()
                   /((2*maxAlphaCo)*mesh.V())
                );
            }
            
            rDeltaT.correctBoundaryConditions();
            
            if (verbose)
            {
                Info<< "Flow time scale min/max = "
                    << gMin(1/rDeltaT.primitiveField())
                    << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
            }
        }
        else
        {
            // Standard Courant number calculation
            #include "CourantNo.H"
            
            // Alpha Courant number calculation
            scalar alphaCoNum = 0.0;
            scalar meanAlphaCoNum = 0.0;
            
            if (mesh.nInternalFaces())
            {
                scalarField sumPhi
                (
                    fvc::surfaceSum(mag(phi))().primitiveField()
                );
                
                alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
                meanAlphaCoNum = 0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
            }
            
            if (verbose)
            {
                Info<< "Interface Courant Number mean: " << meanAlphaCoNum
                    << " max: " << alphaCoNum << endl;
            }
            
            // Standard time step setting
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update laser model
        laser.update();
        laser.correct(runTime.value());
        tmp<volScalarField> laserSrc = laser.source();
        mixture.setQLaser(laserSrc());
        // Refresh cached two-temperature properties in case the model updated
        mixture.setClTTM(ttm.Cl());

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            bool alphaSubCycleExecuted = false;
           bool interfaceCorrectionAppliedInSubCycle = false;

            // Update mixture properties and phase-change sources for the
            // alpha equation
            mixture.correct();
            const volScalarField& phaseChangeSource = mixture.phaseChangeSource();
            const volScalarField& phaseChangeRelaxCoeff =
                mixture.phaseChangeRelaxCoeff();            
            const volScalarField& dgdt = mixture.dgdt();

            // Solve alpha transport using the unified compressible path
#include "compressibleAlphaEqnSubCycle.H"

transportModel.correctPhasePhi();
mixture.correct();

            // Only complain while the laser window is active
            const scalar tnow = runTime.value();
            if (tnow >= laser.laserStartTime() && tnow <= laser.laserEndTime())
            {
                const dimensionedScalar maxQL = max(laserSrc());   // Q_laser field
                const scalar eps = 1e-6;

                if (maxQL.value() < eps)
                {
                    WarningInFunction
                        << "Maximum Q_laser (" << maxQL.value()
                        << ") below threshold " << eps << endl;
                }

                if (verbose)
                {
                    Info<< "max(Q_laser) = " << maxQL.value()
                        << ", max(Tl_) = " << max(ttm.Tl()).value() << nl;
                }
            }
            ttm.solve(laserSrc(), phaseChangeSource, phaseChangeRelaxCoeff);

#include "TEqn.H"

// Recoil update: only if alpha subcycle didn’t run
if (!alphaSubCycleExecuted && pInterfaceCapturing.valid())
{
    const label interval = pInterfaceCapturing->recoilUpdateInterval();
    if (interval <= 1 || recoilCallCount++ % interval == 0)
    {
        pInterfaceCapturing->calculateRecoilPressure();
    }
}

// Apply interface-capturing corrections after recoil update
if
(
    useAdvancedCapturing
 && pInterfaceCapturing.valid()
 && !interfaceCorrectionAppliedInSubCycle
)
{
    pInterfaceCapturing->correct();
}

if (verbose)
{
    Info<< "max recoilPressure = "
        << max(pInterfaceCapturing->recoilPressure()).value() << " Pa" << endl;
}

#include "UEqn.H"
            
            
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
            
            if (pimple.turbCorr())
            {
                transportModel.correct();
            }
        }

       // Compute domain-integrated energy components
        dimensionedScalar Ek = fvc::domainIntegrate
        (
            0.5*rho*magSqr(U)
        );

        const dimensionedScalar& Ce_ = ttm.Ce();
        const dimensionedScalar& Cl_ = ttm.Cl();

        dimensionedScalar Ee = fvc::domainIntegrate
        (
            Ce_*ttm.Te()
        );

        dimensionedScalar Elattice = fvc::domainIntegrate
        (
            Cl_*ttm.Tl()
        );

        const dimensionedScalar L = mixture.latentHeat();
        dimensionedScalar El = fvc::domainIntegrate
        (
            rho1*L*alpha1
        );

        dimensionedScalar Etot = Ek + Ee + Elattice + El;

        static scalar prevEtot = Etot.value();
        scalar dE = Etot.value() - prevEtot;
        scalar relChange = mag(Etot.value() - prevEtot)/max(mag(prevEtot), VSMALL);
        
        if (verbose)
        {
            Info<< "Total energy change: " << dE << " J (" << relChange << ")" << endl;
        }

        if (relChange > energyTolerance)
        {
            WarningInFunction
                << "    Relative energy change " << relChange
                << "    exceeds energyTolerance (" << energyTolerance << ")" << nl
                << "    Ek = " << Ek.value() << " J" << nl
                << "    Ee = " << Ee.value() << " J" << nl
                << "    Elattice = " << Elattice.value() << " J" << nl
                << "    Elatent = " << El.value() << " J" << nl
                << "    Etot = " << Etot.value() << " J" << endl;
        }
                prevEtot = Etot.value();

        // Write additional model fields and data
        if (runTime.writeTime())
        {
            laser.write();
            ttm.write();
            if (pInterfaceCapturing.valid())
            {
                pInterfaceCapturing->write();
            }
        }

        runTime.write();
        if (!U.headerOk())
        {
            WarningInFunction
                << "Field U failed to write. Ensure initial conditions"
                << " and write intervals are set." << endl;
        }

        if (!p_rgh.headerOk())
        {
            WarningInFunction
                << "Field p_rgh failed to write. Ensure initial conditions"
                << " and write intervals are set." << endl;
        }

        runTime.printExecutionTime(Info);
    }
    
    // Clean up dynamically allocated interface capturing object
    if (pInterfaceCapturing.valid())
    {
        if (verbose)
        {
            Info<< "Interface capturing object will be automatically cleaned up" << endl;
        }
        pInterfaceCapturing->clearInternalFields();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
