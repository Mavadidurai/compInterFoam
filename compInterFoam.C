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
#include "dynamicFvMesh.H"
#include "OSspecific.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "autoPtr.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "fvcSmooth.H"
#include "twoPhaseMixtureThermo.H"
#include "extrapolatedCalculatedFvPatchFields.H"
// LIFT-specific model headers
#include "femtosecondLaserModel.H"
#include "twoTemperatureModel.H"
#include "advancedInterfaceCapturing.H"
#include "enhancedLIFTPhysics.H"
#include "OFstream.H"
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>

Foam::Switch verbose(false);
namespace
{
    inline void refreshLegacyRecoilPressure(Foam::volScalarField* fieldPtr)
    {
        if (!fieldPtr)
        {
            return;
        }

        Foam::volScalarField& recoilField = *fieldPtr;
        recoilField.primitiveFieldRef() = 0.0;
        recoilField.correctBoundaryConditions();
    }
    inline Foam::dictionary compInterFoamCoeffsDict(const Foam::fvMesh& mesh)
    {
        const Foam::dictionary& solutionDict = mesh.solutionDict();

        if (solutionDict.found("compInterFoamCoeffs"))
        {
            return solutionDict.subDict("compInterFoamCoeffs");
        }

        return Foam::dictionary();
    }
    static void liftProcessTracker
    (
        const Foam::Time& runTime,
        const Foam::fvMesh& mesh,
        const Foam::twoTemperatureModel& ttm,
        const Foam::volScalarField& alpha1,
        const Foam::volScalarField& rho,
        const Foam::volScalarField& p_rgh,
        const Foam::volScalarField& gh,
        const Foam::twoPhaseMixtureThermo& mixture,
        const Foam::volVectorField& U,
        const Foam::volScalarField* recoilPressurePtr,
        const Foam::Switch trackerEnabled
    )
    {
        static bool statusPrinted = false;

        if (!trackerEnabled)
        {
            if (!statusPrinted && Foam::Pstream::master())
            {
                Info<< "Lift process tracker disabled by controlDict entry"
                    << " (enableLiftProcessTracker = false)" << Foam::endl;
                statusPrinted = true;
            }
            return;
        }

        if (!statusPrinted && Foam::Pstream::master())
        {
            Info<< "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << Foam::nl
                << "   ADVANCED LIFT PROCESS TRACKER v2.0" << Foam::nl
                << "   Physics-based diagnostics with real-time prediction" << Foam::nl
                << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << Foam::endl;
            statusPrinted = true;
        }

        const Foam::scalar t = runTime.value();
        const Foam::volScalarField& Te = ttm.Te();
        const Foam::volScalarField& Tl = ttm.Tl();
        const Foam::volScalarField& QLaser = mixture.Q_laser();

        Foam::tmp<Foam::volScalarField> pTmp(p_rgh + rho*gh);
        const Foam::volScalarField& p = pTmp();

        const Foam::scalar maxTe = Foam::gMax(Te);
        const Foam::scalar maxTl = Foam::gMax(Tl);
        const Foam::scalar minTl = Foam::gMin(Tl);
        const Foam::scalar avgTl = Foam::gAverage(Tl);
        const Foam::scalar maxRecoil = recoilPressurePtr ? Foam::gMax(*recoilPressurePtr) : 0.0;
        const Foam::scalar minRecoil = recoilPressurePtr ? Foam::gMin(*recoilPressurePtr) : 0.0;
        const Foam::scalar peakQLaser = Foam::gMax(QLaser);
        const Foam::scalar laserPowerMetal =
            fvc::domainIntegrate(alpha1*QLaser).value();
        const Foam::scalar maxPressure = Foam::gMax(p);
        const Foam::scalar minPressure = Foam::gMin(p);
        const Foam::scalar avgPressure = Foam::gAverage(p);
        const Foam::scalar tempSpread = maxTl - minTl;

        auto scalarValue = [&](Foam::scalar val, const std::string& unit)
        {
            std::ostringstream valueStream;
            const Foam::scalar absVal = Foam::mag(val);

            if (absVal >= 1e4 || (absVal > Foam::SMALL && absVal <= 1e-3))
            {
                valueStream << std::scientific << std::setprecision(2) << val;
            }
            else
            {
                unsigned int precision = 2;

                if (absVal < 1.0)
                {
                    precision = (absVal < 0.01) ? 4 : 3;
                }

                valueStream << std::fixed << std::setprecision(precision) << val;
            }

            if (!unit.empty())
            {
                valueStream << ' ' << unit;
            }

            return valueStream.str();
        };
        
        const Foam::scalarField& cellVolumes = mesh.V();
        Foam::scalarField alphaVol(alpha1.primitiveField());
        forAll(alphaVol, cellI)
        {
            alphaVol[cellI] *= cellVolumes[cellI];
        }
        Foam::scalar metalVol = Foam::gSum(alphaVol)*1e18;

        static Foam::scalar initialMetalVol = -1.0;
        static Foam::scalar peakMetalVol = -1.0;
        static bool missingStartFieldWarned = false;

        if (initialMetalVol < 0)
        {
            Foam::IOobject alphaInitialIO
            (
                alpha1.name(),
                runTime.timeName(runTime.startTime().value()),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            );

            if (alphaInitialIO.typeHeaderOk<Foam::volScalarField>())
            {
                Foam::volScalarField alphaInitial(alphaInitialIO, mesh);
                Foam::scalarField alphaInitialVol(alphaInitial.primitiveField());
                forAll(alphaInitialVol, cellI)
                {
                    alphaInitialVol[cellI] *= cellVolumes[cellI];
                }

                initialMetalVol = Foam::gSum(alphaInitialVol)*1e18;
                peakMetalVol = initialMetalVol;
            }
            else
            {
                // Regression: continue tracking when the start-time directory
                // was purged and the initial field is unavailable.
                initialMetalVol = metalVol;
                peakMetalVol = metalVol;

                if (!missingStartFieldWarned)
                {
                    missingStartFieldWarned = true;

                    if (Foam::Pstream::master())
                    {
                        WarningInFunction
                            << "Start-time field '" << alpha1.name() << "' not found in "
                            << runTime.timeName(runTime.startTime().value())
                            << "; using current field as initial state." << Foam::endl;
                    }
                }
            }
        }

        if (metalVol > peakMetalVol)
        {
            peakMetalVol = metalVol;
        }

        Foam::scalar metalLoss = 0.0;

        if (initialMetalVol > Foam::SMALL)
        {
            metalLoss = (initialMetalVol - metalVol)/initialMetalVol*100.0;
            metalLoss = Foam::max(metalLoss, Foam::scalar(0));
        }

        const Foam::scalarField& alphaInternal = alpha1.primitiveField();
        const Foam::vectorField& UInternal = U.primitiveField();

        Foam::scalar maxVel = 0.0;
        Foam::scalar weightedVelSum = 0.0;
        Foam::scalar metalVolumeWeight = 0.0;

        forAll(alphaInternal, cellI)
        {
            if (alphaInternal[cellI] > 0.01)
            {
                const Foam::scalar velMag = Foam::mag(UInternal[cellI]);
                maxVel = Foam::max(maxVel, velMag);
                const Foam::scalar weight = alphaInternal[cellI]*cellVolumes[cellI];
                weightedVelSum += velMag*weight;
                metalVolumeWeight += weight;
            }
        }

        Foam::reduce(maxVel, Foam::maxOp<Foam::scalar>());
        Foam::reduce(weightedVelSum, Foam::sumOp<Foam::scalar>());
        Foam::reduce(metalVolumeWeight, Foam::sumOp<Foam::scalar>());

        Foam::scalar avgVel = 0.0;
        if (metalVolumeWeight > Foam::VSMALL)
        {
            avgVel = weightedVelSum/metalVolumeWeight;
        }

        // ========== ADVANCED DIAGNOSTICS ==========

        // 1. Vapor plume tracking
        Foam::scalar vaporVolume = 0.0;
        Foam::scalar vaporMaxVel = 0.0;
        Foam::scalar vaporAvgVel = 0.0;
        Foam::scalar vaporWeightSum = 0.0;
        Foam::label vaporCells = 0;

        const Foam::scalar alpha2Threshold = 0.99; // Almost pure gas
        const Foam::scalar rhoGasThreshold = 10.0; // High-density gas (plume)

        forAll(alphaInternal, cellI)
        {
            if (alphaInternal[cellI] < (1.0 - alpha2Threshold)) // Gas phase
            {
                const Foam::scalar rhoCell = rho.primitiveField()[cellI];
                const Foam::scalar TCell = Tl.primitiveField()[cellI];

                // Detect high-density/high-temp gas (vapor plume, not ambient)
                if (rhoCell > rhoGasThreshold || TCell > 500.0)
                {
                    vaporCells++;
                    const Foam::scalar cellVol = cellVolumes[cellI];
                    vaporVolume += cellVol;

                    const Foam::scalar velMag = Foam::mag(UInternal[cellI]);
                    vaporMaxVel = Foam::max(vaporMaxVel, velMag);

                    const Foam::scalar weight = rhoCell * cellVol;
                    vaporWeightSum += velMag * weight;
                    vaporAvgVel += weight;
                }
            }
        }

        Foam::reduce(vaporVolume, Foam::sumOp<Foam::scalar>());
        Foam::reduce(vaporMaxVel, Foam::maxOp<Foam::scalar>());
        Foam::reduce(vaporWeightSum, Foam::sumOp<Foam::scalar>());
        Foam::reduce(vaporAvgVel, Foam::sumOp<Foam::scalar>());
        Foam::reduce(vaporCells, Foam::sumOp<Foam::label>());

        if (vaporAvgVel > Foam::VSMALL)
        {
            vaporAvgVel = vaporWeightSum / vaporAvgVel;
        }

        // 2. Interface topology analysis
        Foam::label interfaceCells = 0;
        Foam::scalar interfaceArea = 0.0;
        const Foam::scalar alphaInterfaceMin = 0.01;
        const Foam::scalar alphaInterfaceMax = 0.99;

        forAll(alphaInternal, cellI)
        {
            const Foam::scalar alpha = alphaInternal[cellI];
            if (alpha > alphaInterfaceMin && alpha < alphaInterfaceMax)
            {
                interfaceCells++;
                // Approximate interface area contribution
                interfaceArea += Foam::pow(cellVolumes[cellI], 2.0/3.0);
            }
        }

        Foam::reduce(interfaceCells, Foam::sumOp<Foam::label>());
        Foam::reduce(interfaceArea, Foam::sumOp<Foam::scalar>());

        // 3. Film separation detection
        static Foam::scalar prevMetalVol = metalVol;
        static Foam::scalar prevAvgVel = avgVel;
        static Foam::scalar prevRecoil = maxRecoil;
        static Foam::scalar maxVelSeen = 0.0;
        static Foam::scalar maxRecoilSeen = 0.0;
        static bool separationDetected = false;

        const Foam::scalar dt = runTime.deltaTValue();
        const Foam::scalar metalLossRate = dt > Foam::SMALL ?
            (prevMetalVol - metalVol) / dt : 0.0;
        const Foam::scalar velAccel = dt > Foam::SMALL ?
            (avgVel - prevAvgVel) / dt : 0.0;

        maxVelSeen = Foam::max(maxVelSeen, avgVel);
        maxRecoilSeen = Foam::max(maxRecoilSeen, maxRecoil);

        // Detect separation: recoil dropping + velocity still high + material loss increasing
        if (!separationDetected && maxRecoilSeen > 1e7 && maxRecoil < 0.5 * maxRecoilSeen && avgVel > 50.0)
        {
            separationDetected = true;

            if (Foam::Pstream::master())
            {
                Info<< "\n"
                    << "╔═══════════════════════════════════════════════════════════════════╗\n"
                    << "║                                                                   ║\n"
                    << "║  🚀🚀🚀  FILM SEPARATION EVENT DETECTED!  🚀🚀🚀                  ║\n"
                    << "║                                                                   ║\n"
                    << "║  Time: " << t*1e12 << " ps                                             ║\n"
                    << "║  Film velocity: " << avgVel << " m/s                                     ║\n"
                    << "║  Recoil pressure drop: " << (1.0 - maxRecoil/maxRecoilSeen)*100.0 << " %                              ║\n"
                    << "║                                                                   ║\n"
                    << "║  Film has detached from donor substrate!                         ║\n"
                    << "║  Now tracking ballistic transfer phase...                        ║\n"
                    << "║                                                                   ║\n"
                    << "╚═══════════════════════════════════════════════════════════════════╝\n"
                    << Foam::endl;
            }
        }

        // 4. Velocity field decomposition
        Foam::scalar velocityStdDev = 0.0;
        Foam::scalar turbulentKE = 0.0;

        if (metalVolumeWeight > Foam::VSMALL)
        {
            forAll(alphaInternal, cellI)
            {
                if (alphaInternal[cellI] > 0.01)
                {
                    const Foam::scalar velMag = Foam::mag(UInternal[cellI]);
                    const Foam::scalar weight = alphaInternal[cellI] * cellVolumes[cellI];
                    const Foam::scalar deviation = velMag - avgVel;
                    velocityStdDev += weight * deviation * deviation;

                    // Turbulent KE ≈ 0.5 * rho * (U - U_mean)^2
                    const Foam::vector velFluctuation = UInternal[cellI] -
                        (UInternal[cellI] / Foam::mag(UInternal[cellI] + Foam::vector::one*Foam::SMALL)) * avgVel;
                    turbulentKE += 0.5 * rho.primitiveField()[cellI] *
                        alphaInternal[cellI] * cellVolumes[cellI] *
                        Foam::magSqr(velFluctuation);
                }
            }

            Foam::reduce(velocityStdDev, Foam::sumOp<Foam::scalar>());
            Foam::reduce(turbulentKE, Foam::sumOp<Foam::scalar>());
            velocityStdDev = Foam::sqrt(velocityStdDev / metalVolumeWeight);
        }

        // 5. Time-to-ejection estimation
        Foam::scalar timeToEjection = -1.0;
        Foam::scalar predictedEjectionVel = -1.0;

        const Foam::scalar ejectionVelThreshold = 100.0; // m/s bulk velocity
        const Foam::scalar criticalRecoil = 1e7; // 10 MPa minimum for ejection

        if (!separationDetected && maxRecoil > criticalRecoil && velAccel > 1e10)
        {
            // Simple linear extrapolation
            const Foam::scalar velToGo = ejectionVelThreshold - avgVel;
            if (velToGo > 0 && velAccel > 0)
            {
                timeToEjection = velToGo / velAccel;
                predictedEjectionVel = avgVel + velAccel * timeToEjection;
            }
        }

        // 6. Energy partitioning
        const Foam::scalar kineticEnergyMetal = 0.5 * metalVolumeWeight * avgVel * avgVel;
        const Foam::scalar totalLaserEnergy = fvc::domainIntegrate(alpha1 * QLaser).value() * t;
        const Foam::scalar kineticEfficiency = totalLaserEnergy > Foam::SMALL ?
            (kineticEnergyMetal / totalLaserEnergy) * 100.0 : 0.0;

        // 7. Shock wave detection
        Foam::scalar maxPressureGradient = 0.0;
        const Foam::vectorField gradP = fvc::grad(p_rgh + rho*gh)().primitiveField();

        forAll(gradP, cellI)
        {
            maxPressureGradient = Foam::max(maxPressureGradient, Foam::mag(gradP[cellI]));
        }
        Foam::reduce(maxPressureGradient, Foam::maxOp<Foam::scalar>());

        const bool shockWavePresent = maxPressureGradient > 1e15; // Pa/m

        // 8. Event notifications
        static bool ejectionPhaseReported = false;
        static bool highVelocityReported = false;
        static bool shockReported = false;

        if (!ejectionPhaseReported && avgVel > ejectionVelThreshold && Foam::Pstream::master())
        {
            ejectionPhaseReported = true;
            Info<< "\n⚡ MILESTONE: Film velocity exceeds ejection threshold ("
                << ejectionVelThreshold << " m/s) at t=" << t*1e12 << " ps\n" << Foam::endl;
        }

        if (!highVelocityReported && avgVel > 500.0 && Foam::Pstream::master())
        {
            highVelocityReported = true;
            Info<< "\n⚡ MILESTONE: High-velocity regime reached (>500 m/s) at t="
                << t*1e12 << " ps\n" << Foam::endl;
        }

        if (!shockReported && shockWavePresent && Foam::Pstream::master())
        {
            shockReported = true;
            Info<< "\n⚡ ALERT: Shock wave detected (∇p = "
                << maxPressureGradient/1e9 << " GPa/m) at t=" << t*1e12 << " ps\n" << Foam::endl;
        }

        prevMetalVol = metalVol;
        prevAvgVel = avgVel;
        prevRecoil = maxRecoil;

        // ========== END ADVANCED DIAGNOSTICS ==========

        if (Foam::Pstream::master())
        {
            Info<< "══════ LIFT STATE SNAPSHOT ══════" << Foam::nl
                << "Time: " << t*1e12 << " ps" << Foam::nl
                << "Temperatures [K]:" << Foam::nl
                << "  max(Te): " << maxTe << Foam::nl
                << "  max(Tl): " << maxTl << Foam::nl
                << "  min(Tl): " << minTl << Foam::nl
                << "  avg(Tl): " << avgTl << Foam::nl
                << "Pressure [MPa]:" << Foam::nl
                << "  min(p): " << minPressure/1e6 << Foam::nl
                << "  avg(p): " << avgPressure/1e6 << Foam::nl
                << "  max(p): " << maxPressure/1e6 << Foam::nl
                << "Velocity [m/s]:" << Foam::nl
                << "  max(|U|): " << maxVel << Foam::nl
                << "  avg(|U|) (metal): " << avgVel << Foam::nl
                << "Recoil pressure [MPa]:" << Foam::nl
                << "  min: " << minRecoil/1e6 << Foam::nl
                << "  max: " << maxRecoil/1e6 << Foam::nl
                << "Laser coupling:" << Foam::nl
                << "  Peak volumetric source: " << peakQLaser << " W/m^3" << Foam::nl
                << "  Metal-absorbed power: " << laserPowerMetal << " W" << Foam::nl
                << "Metal volume: " << metalVol << " µm³" << Foam::nl
                << "Metal loss: " << metalLoss << " %" << Foam::nl
                << "════════════════════════════════" << Foam::endl;
        }
        static Foam::label currentPhase = 0;
        static bool phaseReported[17] = {};
        static bool initialised = false;
        static bool heatingObserved = false;
        static bool latticeExceededAmbient = false;

        const bool laserActive = (peakQLaser > 1e6);
        const bool electronHeated = (maxTe > 500);
        const bool latticeHeating = (maxTl > 500);
        const bool approachingMelt = (maxTl > 1500);
        const bool melting = (maxTl > 1941);
        const bool vaporizing = (maxTl > 3560);
        const bool plasmaConditions = (maxTl > 5000 && maxPressure > 1e7);
        const bool recoilSignificant = (maxRecoil > 1e7);
        const bool ejecting = (maxVel > 10);
        const bool accelerating = (maxVel > 100);
        const bool significantTransfer = (metalLoss > 2.0);

        if (!heatingObserved && (electronHeated || maxTl > 320))
        {
            heatingObserved = true;
        }

        if (!latticeExceededAmbient && maxTl > 500)
        {
            latticeExceededAmbient = true;
        }

        const bool coolingStarted = (maxTl < 3000 && t > 5e-11);
        const bool belowVapor = (maxTl < 3560);
        const bool belowMelt = (maxTl < 1941);
        const bool nearAmbient = (heatingObserved && latticeExceededAmbient && coolingStarted && maxTl < 500);
        const bool thermallyStable = (nearAmbient && tempSpread < 200);

        Foam::label detectedPhase = 0;

        if (laserActive) detectedPhase = Foam::max(detectedPhase, Foam::label(1));
        if (electronHeated) detectedPhase = Foam::max(detectedPhase, Foam::label(2));
        if (latticeHeating) detectedPhase = Foam::max(detectedPhase, Foam::label(3));
        if (approachingMelt) detectedPhase = Foam::max(detectedPhase, Foam::label(4));
        if (melting) detectedPhase = Foam::max(detectedPhase, Foam::label(5));
        if (vaporizing) detectedPhase = Foam::max(detectedPhase, Foam::label(6));
        if (plasmaConditions) detectedPhase = Foam::max(detectedPhase, Foam::label(7));
        if (recoilSignificant) detectedPhase = Foam::max(detectedPhase, Foam::label(8));
        if (ejecting) detectedPhase = Foam::max(detectedPhase, Foam::label(9));
        if (accelerating) detectedPhase = Foam::max(detectedPhase, Foam::label(10));
        if (significantTransfer) detectedPhase = Foam::max(detectedPhase, Foam::label(11));
        if (coolingStarted) detectedPhase = Foam::max(detectedPhase, Foam::label(12));
        if (belowVapor && coolingStarted) detectedPhase = Foam::max(detectedPhase, Foam::label(13));
        if (belowMelt && coolingStarted) detectedPhase = Foam::max(detectedPhase, Foam::label(14));
        if (nearAmbient) detectedPhase = Foam::max(detectedPhase, Foam::label(15));
        if (thermallyStable && nearAmbient) detectedPhase = Foam::max(detectedPhase, Foam::label(16));

        detectedPhase = Foam::max(detectedPhase, currentPhase);

        if (!initialised)
        {
            currentPhase = detectedPhase;
            for (Foam::label i = 0; i <= currentPhase && i < 17; ++i)
            {
                phaseReported[i] = true;
            }
            initialised = true;
        }

        if (detectedPhase > currentPhase)
        {
            for (Foam::label phase = currentPhase + 1; phase <= detectedPhase && phase < 17; ++phase)
            {
                currentPhase = phase;
                if (phaseReported[currentPhase])
                {
                    continue;
                }

                phaseReported[currentPhase] = true;

                if (Foam::Pstream::master())
                {
                    Info<< '\n' << std::string(70, '=') << Foam::nl;

                    switch (currentPhase)
                    {
                        case 1:
                            Info<< ">>> PHASE 1: LASER ACTIVATION" << Foam::nl
                                << "    Laser power: " << scalarValue(laserPowerMetal, "W") << Foam::nl
                                << "    Peak volumetric source: " << scalarValue(peakQLaser, "W/m³") << Foam::endl;
                            break;
                        case 2:
                            Info<< ">>> PHASE 2: ELECTRON ABSORPTION" << Foam::nl
                                << "    Electron temp: " << maxTe << " K" << Foam::nl
                                << "    Lattice temp: " << maxTl << " K" << Foam::endl;
                            break;
                        case 3:
                            Info<< ">>> PHASE 3: LATTICE HEATING" << Foam::nl
                                << "    Lattice temp: " << maxTl << " K" << Foam::nl
                                << "    Heating rate active" << Foam::endl;
                            break;
                        case 4:
                            Info<< ">>> PHASE 4: APPROACHING MELT" << Foam::nl
                                << "    Lattice temp: " << maxTl << " K" << Foam::nl
                                << "    Near melting point (1941K)" << Foam::endl;
                            break;
                        case 5:
                            Info<< ">>> PHASE 5: MELTING" << Foam::nl
                                << "    Lattice temp: " << maxTl << " K" << Foam::nl
                                << "    Phase change: solid → liquid" << Foam::endl;
                            break;
                        case 6:
                            Info<< ">>> PHASE 6: VAPORIZATION" << Foam::nl
                                << "    Lattice temp: " << maxTl << " K" << Foam::nl
                                << "    Phase change: liquid → vapor" << Foam::endl;
                            break;
                        case 7:
                            Info<< ">>> PHASE 7: PLASMA FORMATION" << Foam::nl
                                << "    Lattice temp: " << maxTl << " K" << Foam::nl
                                << "    Pressure: " << maxPressure/1e6 << " MPa" << Foam::nl
                                << "    High-density plasma plaque" << Foam::endl;
                            break;
                        case 8:
                            Info<< ">>> PHASE 8: RECOIL PRESSURE" << Foam::nl
                                << "    Recoil: " << maxRecoil/1e6 << " MPa" << Foam::nl
                                << "    Momentum transfer initiated" << Foam::endl;
                            break;
                        case 9:
                            Info<< ">>> PHASE 9: EJECTION START" << Foam::nl
                                << "    Max velocity: " << maxVel << " m/s" << Foam::nl
                                << "    Metal loss: " << metalLoss << "%" << Foam::endl;
                            break;
                        case 10:
                            Info<< ">>> PHASE 10: ACCELERATION" << Foam::nl
                                << "    Max velocity: " << maxVel << " m/s" << Foam::nl
                                << "    Rapid momentum transfer" << Foam::endl;
                            break;
                        case 11:
                            Info<< ">>> PHASE 11: MATERIAL TRANSFER" << Foam::nl
                                << "    Metal loss: " << metalLoss << "%" << Foam::nl
                                << "    Droplet formation" << Foam::endl;
                            break;
                        case 12:
                            Info<< ">>> PHASE 12: COOLING START" << Foam::nl
                                << "    Max temp: " << maxTl << " K" << Foam::nl
                                << "    Thermal decay initiated" << Foam::endl;
                            break;
                        case 13:
                            Info<< ">>> PHASE 13: BELOW VAPORIZATION" << Foam::nl
                                << "    Max temp: " << maxTl << " K" << Foam::nl
                                << "    Condensation occurring" << Foam::endl;
                            break;
                        case 14:
                            Info<< ">>> PHASE 14: SOLIDIFICATION" << Foam::nl
                                << "    Max temp: " << maxTl << " K" << Foam::nl
                                << "    Phase change: liquid → solid" << Foam::endl;
                            break;
                        case 15:
                            Info<< ">>> PHASE 15: NEAR AMBIENT" << Foam::nl
                                << "    Max temp: " << maxTl << " K" << Foam::nl
                                << "    Approaching equilibrium" << Foam::endl;
                            break;
                        case 16:
                            Info<< ">>> PHASE 16: THERMAL EQUILIBRIUM" << Foam::nl
                                << "    Temp spread: " << tempSpread << " K" << Foam::nl
                                << "    LIFT PROCESS COMPLETE" << Foam::endl;
                            break;
                        default:
                            break;
                    }

                    Info<< "    Time: " << t*1e12 << " ps" << Foam::nl;
                    Info<< std::string(70, '=') << "\n" << Foam::endl;
                }
            }
        }

        static Foam::label reportCounter = 0;
        ++reportCounter;

        if
        (
            (reportCounter == 1 || reportCounter % 10 == 0)
         && Foam::Pstream::master()
        )

        {
            Foam::word phaseName("UNKNOWN");
            switch (currentPhase)
            {
                case 0: phaseName = "IDLE"; break;
                case 1: phaseName = "LASER_ACTIVE"; break;
                case 2: phaseName = "ABSORPTION"; break;
                case 3: phaseName = "HEATING"; break;
                case 4: phaseName = "PRE-MELT"; break;
                case 5: phaseName = "MELTING"; break;
                case 6: phaseName = "VAPORIZING"; break;
                case 7: phaseName = "PLASMA"; break;
                case 8: phaseName = "RECOIL"; break;
                case 9: phaseName = "EJECTION"; break;
                case 10: phaseName = "ACCELERATION"; break;
                case 11: phaseName = "TRANSFER"; break;
                case 12: phaseName = "COOLING"; break;
                case 13: phaseName = "CONDENSING"; break;
                case 14: phaseName = "SOLIDIFYING"; break;
                case 15: phaseName = "NEAR_AMBIENT"; break;
                case 16: phaseName = "EQUILIBRIUM"; break;
            }

            Foam::scalar progressIndex = Foam::scalar(currentPhase);

            if (currentPhase == 3)
            {
                const Foam::scalar heatingStart = 500.0;
                const Foam::scalar heatingTarget = 1500.0;

                if (heatingTarget > heatingStart)
                {
                    Foam::scalar heatingFraction = (maxTl - heatingStart)/(heatingTarget - heatingStart);
                    heatingFraction = Foam::max(Foam::scalar(0), Foam::min(heatingFraction, Foam::scalar(0.999)));
                    progressIndex = Foam::scalar(currentPhase) + heatingFraction;
                }
            }

             const Foam::scalar progress = Foam::min(progressIndex/16.0*100.0, Foam::scalar(100));

            const int frameWidth = 80;
            const int innerWidth = frameWidth - 2;
            const int labelWidth = 30;
            const int valueWidth = innerWidth - labelWidth - 2;

            auto centeredLine = [&](const std::string& text)
            {
                const std::size_t textLen = text.size();
                const int paddingTotal = innerWidth - static_cast<int>(textLen);
                const int paddingLeft = Foam::max(paddingTotal/2, 0);
                const int paddingRight = Foam::max(paddingTotal - paddingLeft, 0);
                Info<< "║" << std::string(paddingLeft, ' ')
                    << text
                    << std::string(paddingRight, ' ') << "║" << Foam::nl;
            };

            auto formattedRow = [&](const std::string& label, const std::string& value)
            {
                std::ostringstream rowStream;
                rowStream << "║ "
                          << std::left << std::setw(labelWidth) << label
                          << std::right << std::setw(valueWidth) << value
                          << " ║";
                Info<< rowStream.str() << Foam::nl;
            };

            // Enhanced progress bar
            auto progressBar = [&](Foam::scalar pct)
            {
                const int barWidth = 50;
                const int filled = static_cast<int>(pct / 100.0 * barWidth);
                std::string bar = "[";
                for (int i = 0; i < barWidth; ++i)
                {
                    if (i < filled) bar += "█";
                    else if (i == filled) bar += "▌";
                    else bar += "░";
                }
                bar += "]";

                // Format progress percentage
                std::ostringstream pctStr;
                pctStr << std::fixed << std::setprecision(1) << pct << "%";

                // Pad to ensure alignment
                while (bar.length() + pctStr.str().length() < 56)
                {
                    bar += " ";
                }
                bar += pctStr.str();

                return bar;
            };

            Info<< "\n╔════════════════════════════════════════════════════════════════════════════════╗" << Foam::nl;
            centeredLine("⚡ ADVANCED LIFT PROCESS STATUS ⚡");
            Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
            formattedRow("Phase:", std::string(phaseName));
            formattedRow("Time:", scalarValue(t*1e12, "ps"));
            Info<< "║ Progress: " << progressBar(progress) << " ║" << Foam::nl;

            // Separation status indicator
            if (separationDetected)
            {
                Info<< "║ " << std::string(76, ' ') << " ║" << Foam::nl;
                centeredLine("🚀 *** FILM SEPARATION DETECTED *** 🚀");
            }

            Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
            centeredLine("🌡️  THERMAL DIAGNOSTICS");
            formattedRow("Electron temp (max):", scalarValue(maxTe, "K"));
            formattedRow("Lattice temp (max):", scalarValue(maxTl, "K"));
            formattedRow("Lattice temp (avg):", scalarValue(avgTl, "K"));
            formattedRow("Temperature spread:", scalarValue(tempSpread, "K"));

            Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
            centeredLine("💨 PRESSURE & DYNAMICS");
            formattedRow("Max pressure:", scalarValue(maxPressure/1e6, "MPa"));
            formattedRow("Recoil pressure:", scalarValue(maxRecoil/1e6, "MPa"));
            formattedRow("Peak recoil seen:", scalarValue(maxRecoilSeen/1e6, "MPa"));
            formattedRow("Pressure gradient:", scalarValue(maxPressureGradient/1e9, "GPa/m"));
            if (shockWavePresent)
            {
                formattedRow("Shock wave:", "DETECTED");
            }

            Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
            centeredLine("🎯 METAL FILM STATUS");
            formattedRow("Volume:", scalarValue(metalVol, "µm³"));
            formattedRow("Volume loss:", scalarValue(metalLoss, "%"));
            formattedRow("Loss rate:", scalarValue(metalLossRate*1e18, "µm³/s"));
            formattedRow("Interface cells:", std::to_string(interfaceCells));
            formattedRow("Interface area:", scalarValue(interfaceArea*1e12, "µm²"));

            Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
            centeredLine("⚡ VELOCITY ANALYSIS");
            formattedRow("Film velocity (avg):", scalarValue(avgVel, "m/s"));
            formattedRow("Film velocity (max):", scalarValue(maxVel, "m/s"));
            formattedRow("Velocity std dev:", scalarValue(velocityStdDev, "m/s"));
            formattedRow("Acceleration:", scalarValue(velAccel/1e12, "km/s²"));
            formattedRow("Turbulent KE:", scalarValue(turbulentKE*1e9, "nJ"));

            Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
            centeredLine("☁️  VAPOR PLUME TRACKING");
            formattedRow("Plume volume:", scalarValue(vaporVolume*1e18, "µm³"));
            formattedRow("Plume cells:", std::to_string(vaporCells));
            formattedRow("Plume velocity (avg):", scalarValue(vaporAvgVel, "m/s"));
            formattedRow("Plume velocity (max):", scalarValue(vaporMaxVel, "m/s"));

            Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
            centeredLine("⚡ ENERGY PARTITIONING");
            formattedRow("Laser power (current):", scalarValue(laserPowerMetal, "W"));
            formattedRow("Peak source:", scalarValue(peakQLaser/1e12, "TW/m³"));
            formattedRow("Metal kinetic energy:", scalarValue(kineticEnergyMetal*1e9, "nJ"));
            formattedRow("Laser→kinetic eff:", scalarValue(kineticEfficiency, "%"));

            if (timeToEjection > 0 && !separationDetected)
            {
                Info<< "╠════════════════════════════════════════════════════════════════════════════════╣" << Foam::nl;
                centeredLine("🔮 EJECTION PREDICTION");
                formattedRow("Est. time to ejection:", scalarValue(timeToEjection*1e12, "ps"));
                formattedRow("Predicted ejection vel:", scalarValue(predictedEjectionVel, "m/s"));
                formattedRow("Time at ejection:", scalarValue((t + timeToEjection)*1e12, "ps"));
            }

            Info<< "╚════════════════════════════════════════════════════════════════════════════════╝" << Foam::nl;

            if (maxTl > 15000)
            {
                Info<< "⚠  WARNING: Temperature exceeds physical limits!" << Foam::endl;
            }
            const Foam::dictionary solverCoeffs = compInterFoamCoeffsDict(mesh);
            Foam::scalar realismLimit = solverCoeffs.lookupOrDefault<Foam::scalar>
            (
                "maxReasonableVelocity",
                0.0
            );

            if (realismLimit <= Foam::SMALL)
            {
                realismLimit = 5000.0;
            }

            if (maxVel > realismLimit)
            {
                Info<< "⚠  WARNING: Velocity exceeds realistic LIFT range!" << Foam::nl
                    << "       Limit: " << realismLimit << " m/s" << Foam::nl
                    << "       Observed: " << maxVel << " m/s" << Foam::nl;

                const Foam::scalar rhoLiquid = mixture.rho1().value();
                Foam::scalar recoilDrivenSpeed = 0.0;

                if (rhoLiquid > Foam::SMALL && maxRecoil > 0)
                {
                    recoilDrivenSpeed = Foam::sqrt(2.0*maxRecoil/rhoLiquid);
                }

                Info<< "       Max recoil pressure: " << maxRecoil/1e6 << " MPa" << Foam::nl;

                if (recoilDrivenSpeed > 0)
                {
                    Info<< "       Ideal recoil jet speed sqrt(2*deltaP/rho): "
                        << recoilDrivenSpeed << " m/s" << Foam::nl;
                }

                Info<< "       Time step: " << runTime.deltaTValue() << " s" << Foam::nl
                    << "       Tl spread: " << tempSpread << " K" << Foam::endl;
            }
            if (metalLoss > 50)
            {
                Info<< "⚠  WARNING: Excessive material loss!" << Foam::endl;
            }
            if (currentPhase < 5 && t > 1e-10)
            {
                Info<< "⚠  WARNING: Process too slow - check laser parameters!" << Foam::endl;
            }

            Info<< Foam::endl;
        }

        static Foam::label csvCounter = 0;
        ++csvCounter;

        if (csvCounter % 25 == 0 && Foam::Pstream::master())
        {
            static bool headerWritten = false;

            Foam::OFstream csvFile
            (
                "liftProcessTracking.csv",
                Foam::IOstreamOption::ASCII,
                Foam::IOstreamOption::UNCOMPRESSED,
                Foam::IOstreamOption::APPEND
            );

            if (!headerWritten)
            {
                csvFile << "time_ps,phase,phase_num,progress_pct,"
                        << "Te_max_K,Tl_max_K,Tl_avg_K,Tl_spread_K,"
                        << "P_max_MPa,recoil_MPa,recoil_max_seen_MPa,"
                        << "metal_vol_um3,metal_loss_pct,metal_loss_rate_um3s,"
                        << "vel_max_ms,vel_avg_ms,vel_stddev_ms,accel_kms2,"
                        << "interface_cells,interface_area_um2,"
                        << "vapor_vol_um3,vapor_cells,vapor_vel_avg_ms,vapor_vel_max_ms,"
                        << "turbulent_KE_nJ,kinetic_eff_pct,"
                        << "pressure_grad_GPam,shock_present,"
                        << "separation_detected,time_to_ejection_ps,"
                        << "laserPower_W,peakQLaser_TWm3" << Foam::endl;
                headerWritten = true;
            }

            Foam::word phaseName("UNKNOWN");
            switch (currentPhase)
            {
                case 0: phaseName = "IDLE"; break;
                case 1: phaseName = "LASER_ACTIVE"; break;
                case 2: phaseName = "ABSORPTION"; break;
                case 3: phaseName = "HEATING"; break;
                case 4: phaseName = "PRE_MELT"; break;
                case 5: phaseName = "MELTING"; break;
                case 6: phaseName = "VAPORIZING"; break;
                case 7: phaseName = "PLASMA"; break;
                case 8: phaseName = "RECOIL"; break;
                case 9: phaseName = "EJECTION"; break;
                case 10: phaseName = "ACCELERATION"; break;
                case 11: phaseName = "TRANSFER"; break;
                case 12: phaseName = "COOLING"; break;
                case 13: phaseName = "CONDENSING"; break;
                case 14: phaseName = "SOLIDIFYING"; break;
                case 15: phaseName = "NEAR_AMBIENT"; break;
                case 16: phaseName = "EQUILIBRIUM"; break;
            }

            Foam::scalar progressIndex = Foam::scalar(currentPhase);

            if (currentPhase == 3)
            {
                const Foam::scalar heatingStart = 500.0;
                const Foam::scalar heatingTarget = 1500.0;

                if (heatingTarget > heatingStart)
                {
                    Foam::scalar heatingFraction = (maxTl - heatingStart)/(heatingTarget - heatingStart);
                    heatingFraction = Foam::max(Foam::scalar(0), Foam::min(heatingFraction, Foam::scalar(0.999)));
                    progressIndex = Foam::scalar(currentPhase) + heatingFraction;
                }
            }

            const Foam::scalar progress = Foam::min(progressIndex/16.0*100.0, Foam::scalar(100));

            csvFile << t*1e12 << ","
                    << phaseName << ","
                    << currentPhase << ","
                    << progress << ","
                    << maxTe << ","
                    << maxTl << ","
                    << avgTl << ","
                    << tempSpread << ","
                    << maxPressure/1e6 << ","
                    << maxRecoil/1e6 << ","
                    << maxRecoilSeen/1e6 << ","
                    << metalVol << ","
                    << metalLoss << ","
                    << metalLossRate*1e18 << ","
                    << maxVel << ","
                    << avgVel << ","
                    << velocityStdDev << ","
                    << velAccel/1e12 << ","
                    << interfaceCells << ","
                    << interfaceArea*1e12 << ","
                    << vaporVolume*1e18 << ","
                    << vaporCells << ","
                    << vaporAvgVel << ","
                    << vaporMaxVel << ","
                    << turbulentKE*1e9 << ","
                    << kineticEfficiency << ","
                    << maxPressureGradient/1e9 << ","
                    << (shockWavePresent ? 1 : 0) << ","
                    << (separationDetected ? 1 : 0) << ","
                    << (timeToEjection > 0 ? timeToEjection*1e12 : -1) << ","
                    << laserPowerMetal << ","
                    << peakQLaser/1e12 << Foam::endl;
        }

        pTmp.clear();
    }
}
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
/*
    const Foam::fileName meshDir = runTime.constant()/"polyMesh";
        if (!Foam::isDir(meshDir) || !Foam::isFile(meshDir/"points"))
    {
        FatalErrorInFunction
            << "Missing initial mesh: " << meshDir/"points" << '\n'
            << "Generate the background mesh (e.g. blockMesh or surface meshing) "
            << "before launching compInterFoam."
            << exit(FatalError);
    }

    autoPtr<dynamicFvMesh> meshPtr(dynamicFvMesh::New(runTime));
    dynamicFvMesh& mesh = meshPtr();
    */
    #include "createTimeControls.H" 
    #ifndef CREATE_FIELDS_DONE
    #include "createFields.H"
    #define CREATE_FIELDS_DONE
    #endif

    const volScalarField& phaseChangeSource = mixture.phaseChangeSource();
    const volScalarField& phaseChangeRelaxCoeff = mixture.phaseChangeRelaxCoeff();
    
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

    const Foam::Switch enableLiftProcessTracker
    (
        runTime.controlDict().lookupOrDefault<Foam::Switch>
        (
            "enableLiftProcessTracker",
            true
        )
    );

    if (Foam::Pstream::master())
    {
        Info<< "Lift process tracker "
            << (enableLiftProcessTracker ? "enabled" : "disabled")
            << " (enableLiftProcessTracker in controlDict)" << Foam::endl;
    }

    runTime.functionObjects().start();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        const bool master = Foam::Pstream::master();
        const polyMesh::readUpdateState meshState = mesh.readUpdate();

        if (meshState != polyMesh::UNCHANGED)
        {
            gh = (g & mesh.C()) - ghRef;
            ghf = (g & mesh.Cf()) - ghRef;
            p_rgh = p - rho*gh;
            p_rgh.correctBoundaryConditions();

            if (meshState == polyMesh::TOPO_CHANGE)
            {
                phi = fvc::flux(U);
                rhoPhi = fvc::interpolate(rho)*phi;
                alphaPhi10 = phi*fvc::interpolate(alpha1);
            }
        }
        #include "readTimeControls.H"

        if (LTS)
        {
            // Traditional inline LTS implementation
            volScalarField& rDeltaT = trDeltaT.ref();
            volScalarField& rSubDeltaT = trSubDeltaT.ref();

            const scalar maxDeltaTValue
            (
                runTime.controlDict().lookupOrDefault<scalar>
                (
                    "maxDeltaT",
                    GREAT
                )
            );
            const scalar rawMinDeltaT =
                runTime.controlDict().lookupOrDefault<scalar>
                (
                    "minDeltaT",
                    0.0
                );
            const scalar minDeltaTValue = rawMinDeltaT > 0 ? rawMinDeltaT : SMALL;
            const dimensionedScalar maxDeltaT
            (
                "maxDeltaT",
                dimTime,
                maxDeltaTValue
            );

            const dimensionedScalar invMaxDeltaT(1/maxDeltaT);

            scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>
            (
                "maxCo",
                0.9
            );
            scalar maxAlphaCo = runTime.controlDict().lookupOrDefault<scalar>
            (
                "maxAlphaCo",
                0.2
            );

            const dictionary& solutionDict = mesh.solutionDict();
            if (solutionDict.found("PIMPLE"))
            {
                const dictionary& pimpleDict = solutionDict.subDict("PIMPLE");
                maxCo = pimpleDict.getOrDefault<scalar>("maxCo", maxCo);
                maxAlphaCo = pimpleDict.getOrDefault<scalar>("maxAlphaCo", maxAlphaCo);
            }

            // Set reciprocal time-step from local Courant number
            rDeltaT.ref() = max
            (
                invMaxDeltaT,
                fvc::surfaceSum(mag(rhoPhi))()
               /((2*maxCo)*mesh.V()*rho)
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
                   *fvc::surfaceSum(mag(phi))()
                   /((2*maxAlphaCo)*mesh.V())
                );
            }

            rDeltaT.correctBoundaryConditions();
            rSubDeltaT = rDeltaT;
            rSubDeltaT = max(invMaxDeltaT, rSubDeltaT);
            rSubDeltaT.correctBoundaryConditions();
            scalar newDeltaT = 1.0/gMax(rDeltaT.primitiveField());
            if (!std::isfinite(newDeltaT) || newDeltaT <= SMALL)
            {
                newDeltaT = Foam::max(runTime.deltaTValue(), minDeltaTValue);
            }

            newDeltaT = Foam::max
            (
                minDeltaTValue,
                Foam::min(maxDeltaTValue, newDeltaT)
            );
            runTime.setDeltaT(newDeltaT);
           
            if (verbose && master)
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
            
            if (verbose && master)
            {
                Info<< "Interface Courant Number mean: " << meanAlphaCoNum
                    << " max: " << alphaCoNum << endl;
            }
            
            // Standard time step setting
            #include "setDeltaT.H"
        }

        ++runTime;
/*        
        static bool dynamicMeshSuppressedPrinted = false;

        if (mesh.dynamic())
        {
            const Switch enableDynamicMeshRefinement =
                runTime.controlDict()
                    .subOrEmptyDict("advancedInterfaceCapturing")
                    .lookupOrDefault<Switch>("enableDynamicMeshRefinement", true);

            if (useAdvancedCapturing && enableDynamicMeshRefinement)
            {
               const bool meshChangedThisStep = mesh.update();

                if (meshChangedThisStep)
                {
                    dynamicMeshSuppressedPrinted = false;

                    if (verbose && master)
                    {
                        Info<< "Dynamic mesh update executed" << endl;
                    }

                    phi = fvc::flux(U);
                    rhoPhi = fvc::interpolate(rho)*phi;
                    alphaPhi10 = phi*fvc::interpolate(alpha1);
                }
            }
            else if (!dynamicMeshSuppressedPrinted && master)
            {
                Info<< "Dynamic mesh refinement suppressed"
                    << (useAdvancedCapturing
                        ? " (enableDynamicMeshRefinement = false in advancedInterfaceCapturing)"
                        : " (useAdvancedInterfaceCapturing = false)")
                    << endl;
                dynamicMeshSuppressedPrinted = true;
            }
        }
*/
        runTime.functionObjects().execute();

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update laser model
        laser.update();
        laser.correct();
        const tmp<volScalarField> tLaserSource = laser.source();
        mixture.setQLaser(tLaserSource());
        const volScalarField& QLaser = mixture.Q_laser();
        if (verbose && master)
        {
            Info<< "debug: max(Q_laser) now = " << gMax(QLaser) << nl;
        }
        // Refresh cached two-temperature properties in case the model updated
        mixture.setClTTM(ttm.Cl());

    // --- Pressure-velocity PIMPLE corrector loop

while (pimple.loop())
{
    mixture.correct();  // Update phase fractions

    #include "compressibleAlphaEqnSubCycle.H"
    transportModel.correctPhasePhi();
    mixture.correct();

    // ✅ UPDATE TEMPERATURES FIRST
    const label nThermalCouplingIter =
        runTime.controlDict().lookupOrDefault<label>("nThermalCouplingIter", 1);

    for (label thermalIter = 0; thermalIter < nThermalCouplingIter; ++thermalIter)
    {
        #include "TEqn.H"

        ttm.solve
        (
            QLaser,
            phaseChangeSource,
            phaseChangeRelaxCoeff,
            gasMetalHeatFlux
        );
        
        mixture.setClTTM(ttm.Cl());
    }

    // Update enhanced LIFT physics (all three models updated in one call)
    if (liftPhysics.valid())
    {
        liftPhysics->updateAll(ttm.Tl(), alpha1, rho, U);

        if (verbose && master && runTime.timeIndex() % 100 == 0)
        {
            Info<< "Enhanced physics updated:" << nl;
            if (liftPhysics->phaseExplosionEnabled())
                Info<< "    Phase explosion active" << nl;
            if (liftPhysics->plasmaEnabled())
                Info<< "    Plasma ionization active" << nl;
            if (liftPhysics->breakupEnabled())
                Info<< "    Droplet breakup active" << nl;
        }
    }
    
    mixture.correct();

    if (useAdvancedCapturing && pInterfaceCapturing.valid())
    {
        pInterfaceCapturing->calculateRecoilPressure();
        
        if (verbose && master && runTime.timeIndex() % 10 == 0)
        {
            Info<< "    Max recoil pressure: "
                << max(pInterfaceCapturing->recoilPressure()).value()/1e6
                << " MPa" << endl;
        }
    }
    else if (legacyRecoilPressure.valid())
    {
        refreshLegacyRecoilPressure(legacyRecoilPressure.ptr());
    }

    // ✅ NOW SOLVE MOMENTUM WITH CORRECT RECOIL
    #include "UEqn.H"

    // Pressure-velocity coupling
    while (pimple.correct())
    {
    if (useAdvancedCapturing && pInterfaceCapturing.valid())
    {
        pInterfaceCapturing->calculateRecoilPressure();
    }
        #include "pEqn.H"
    }

    if (pimple.turbCorr())
    {
        transportModel.correct();
    }
}
// Apply droplet breakup to alpha field (interface remapping)
if (liftPhysics.valid() && liftPhysics->breakupEnabled())
{
    liftPhysics->applyBreakup(alpha1, U);
}
        dimensionedScalar Ek = fvc::domainIntegrate(0.5*rho*magSqr(U));
        tmp<volScalarField> tCe = ttm.electronHeatCapacity();
        const volScalarField& CeField = tCe();
        const volScalarField& TeField = ttm.Te();
        const volScalarField& TlField = ttm.Tl();
        const dimensionedScalar& Cl_ = ttm.Cl();

        dimensionedScalar Ee = fvc::domainIntegrate(alpha1*CeField*TeField);
        dimensionedScalar Elattice = fvc::domainIntegrate(alpha1*Cl_*TlField);
        const dimensionedScalar El = ttm.cumulativePhaseChangeEnergy();
        const volScalarField& he2 = mixture.thermo2().he();
        dimensionedScalar Egas =
            fvc::domainIntegrate(alpha2*rho2*he2 - alpha2*p);

        dimensionedScalar Etot = Ek + Ee + Elattice + El + Egas;

        // Track cumulative boundary energy flux to account for energy leaving domain
        static dimensionedScalar cumulativeBoundaryFlux
        (
            "cumulativeBoundaryFlux",
            dimEnergy,
            0.0
        );
        // Calculate energy flux through boundaries this timestep
        // Energy flux = (kinetic + enthalpy) * mass flux = (0.5*U^2 + he) * rho * phi
        tmp<volScalarField> tSpecificEnergy(0.5*magSqr(U) + he2);
        const volScalarField& specificEnergy = tSpecificEnergy();
        surfaceScalarField energyFlux
        (
            "energyFlux",
            fvc::interpolate(rho*specificEnergy)*phi
        );
        // Integrate boundary fluxes (positive = leaving domain)
        dimensionedScalar boundaryEnergyThisStep
        (
            "boundaryEnergyThisStep",
            dimPower,
            0.0
        );
        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            // Only account for patches where material can leave domain
            // Skip walls, symmetry, empty, and processor boundaries
            const word& patchType = patch.type();
            const bool isFlowBoundary =
                (patchType != "wall" && patchType != "symmetryPlane"
                 && patchType != "empty" && patchType != "processor"
                 && patchType != "cyclic");

            if (isFlowBoundary)
            {
                const scalarField& patchFlux = energyFlux.boundaryField()[patchI];
                boundaryEnergyThisStep.value() +=
                    gSum(patchFlux);
            }
        }
        cumulativeBoundaryFlux += boundaryEnergyThisStep*dimensionedScalar("dt", dimTime, runTime.deltaT().value());
        static scalar prevEtot = Etot.value();
        const scalar minEnergyForCheck =
            runTime.controlDict().lookupOrDefault<scalar>
            (
                "energyCheckMinEnergy",
                1e-6
            );

        const scalar prevEtotMag = mag(prevEtot);
        const scalar currEtotMag = mag(Etot.value());

        const scalar denom = max(max(prevEtotMag, currEtotMag), minEnergyForCheck);
        // Account for boundary losses in energy change calculation

        const scalar energyChangeWithFlux =

            mag(Etot.value() - prevEtot + cumulativeBoundaryFlux.value());

        const scalar relChange = energyChangeWithFlux/denom;
        const bool accountForBoundaryFlux =

            runTime.controlDict().lookupOrDefault<Switch>
            (
                "energyCheckIncludeBoundaryFlux",
                true
            );
        if
        (
            prevEtotMag > minEnergyForCheck
         && relChange > energyTolerance.value()
         && accountForBoundaryFlux
        )
        {
            WarningInFunction
                << "Relative energy change " << relChange
                << " exceeds energyTolerance (" << energyTolerance.value() << ")" << nl
                << "Domain energy: " << Etot.value() << " J" << nl
                << "Previous: " << prevEtot << " J" << nl
                << "Cumulative boundary loss: " << cumulativeBoundaryFlux.value() << " J" << nl
                << "Energy change (with boundary): " << energyChangeWithFlux << " J"
                << endl;
        }
        else if
        (
            prevEtotMag > minEnergyForCheck
         && mag(Etot.value() - prevEtot)/denom > energyTolerance.value()
         && !accountForBoundaryFlux
        )
        {
            WarningInFunction
                << "Relative energy change " << mag(Etot.value() - prevEtot)/denom
                << " exceeds energyTolerance (" << energyTolerance.value() << ")" << nl
                << "(Boundary flux accounting disabled)" << endl;
        }
            if (verbose && master)
        {
            Info<< "Energy totals [J]: Ek=" << Ek.value()
                << " Ee=" << Ee.value()
                << " Elattice=" << Elattice.value()
                << " Elatent=" << El.value()
                << " Egas=" << Egas.value()
                << " Etot=" << Etot.value() << nl
                << " Boundary flux [J]: This step=" << boundaryEnergyThisStep.value()
                << " Cumulative=" << cumulativeBoundaryFlux.value() << endl;
        }

        prevEtot = Etot.value();
        tCe.clear();

        liftProcessTracker
        (
            runTime,
            mesh,
            ttm,
            alpha1,
            rho,
            p_rgh,
            gh,
            mixture,
            U,
            recoilPressurePtr,
            enableLiftProcessTracker
        );

        runTime.functionObjects().execute();

        // Write additional model fields and data
        if (runTime.write())
        {
            laser.write();
            ttm.write();
            if (useAdvancedCapturing && pInterfaceCapturing.valid())
            {
                pInterfaceCapturing->write();
            }
            else if (legacyRecoilPressure.valid())
            {
                legacyRecoilPressure->write();
            }
        }

        runTime.printExecutionTime(Info);
    }
    
    // Clean up dynamically allocated interface capturing object
    if (useAdvancedCapturing && pInterfaceCapturing.valid())
    {
        if (verbose && Foam::Pstream::master())
        {
            Info<< "Interface capturing object will be automatically cleaned up" << endl;
        }
        pInterfaceCapturing->clearInternalFields();
    }

    Info<< "End\n" << endl;

    return 0;
}
