/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024
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

Class
    Foam::advancedInterfaceCapturing

Description
    Helper class for applying recoil pressure limits and diagnostics for
    Laser-Induced Forward Transfer (LIFT) simulations. Supports optional
    ramping and per-step change limiting of the recoil load to avoid sudden
    jumps during start-up.
\*---------------------------------------------------------------------------*/
#ifndef advancedInterfaceCapturing_H
#define advancedInterfaceCapturing_H
#include "fvMesh.H"
#include "volFields.H"
#include "fvc.H"
#include "twoPhaseMixtureThermo.H"
namespace Foam
{
class advancedInterfaceCapturing
{
private:
    // References
    const fvMesh& mesh_;
    volScalarField& alpha1_;         // Phase fraction
    const twoPhaseMixtureThermo& mixture_;
    const volScalarField& T_;        // Gas temperature field (fallback)
    const volScalarField* TlPtr_;    // Optional lattice temperature field (falls back to T_)
    // Phase-change parameters provided by twoPhaseMixtureThermo
    dimensionedScalar meltingTemp_;             // Melting temperature (Tsol)
    dimensionedScalar vaporTemp_;               // Vaporization temperature (Tvap)
        // User-configurable recoil pressure parameters
    dimensionedScalar pressureScale_; // Scaling factor for evaporation rate [Pa·s]
    scalar recoilMax_;               // Base recoil pressure cap
    dimensionedScalar phaseChangeTempOffset_;   // Offset for phase-change threshold
    dimensionedScalar maxPhysicalTemperature_;  // Optional physical temperature clamp
    Switch clampRecoil_;             // Enable/disable recoil pressure clamping
    Switch scaleRecoilMax_;          // Scale recoilMax with interface weighting
    scalar recoilRelax_;             // Temporal smoothing factor for recoil pressure
    // Valid alpha range (controlDict → advancedInterfaceCapturing)
    scalar alphaMin_;               // Minimum alpha for recoil activation
    scalar alphaMax_;               // Maximum alpha for recoil activation
    scalar massRateEps_;            // Tolerance for evaporation rate comparisons
    scalar metalAlphaCutoff_;       // Minimum alpha treated as metal-dominated
    // Field for recoil pressure
    volScalarField recoilPressure_;  // Recoil pressure from vaporization
    scalarField previousRecoilPressure_;        // Cached recoil pressure for smoothing
    bool havePreviousRecoil_;                    // Flag indicating cache validity
    scalar rampProgress_;            // Progress of recoil ramp (0-1)
    scalar rampIncrement_;           // Increment per time step for recoil ramp
    scalar recoilMaxDelta_;          // Maximum recoil change rate [Pa/s]
    scalar recoilSmoothCoeff_;       // Weight for optional smoothing
    label  recoilSmoothIters_;       // Number of smoothing iterations
    bool femtosecondCase_;           // Flag to relax damping for femtosecond runs
    dimensionedScalar boltzmannConstant_;      // Boltzmann constant
    dimensionedScalar vaporParticleMass_;      // Mass of vapor particle
    scalar momentumAccommodationCoeff_;        // Momentum accommodation coefficient
    scalar recoilCoefficient_;       // p_recoil/p_sat ratio (Anisimov ~0.54)
    label relaxTimeIndex_;           // Time index of stored relaxation reference
    // Confined vapour pocket (approach A): 0-D state, global (reduced)
    scalar pocketMass_;              // evaporated mass in the pocket [kg] (wedge)
    scalar pocketDisp_;              // mean displacement of the film away from the glass [m]
    scalar pocketRadius_;            // largest radius that has evaporated on the glass side [m]
    scalar pocketPressure_;          // current pocket pressure [Pa]
    label  pocketTimeIndex_;         // time index of the last state update
    scalar pocketTime_;              // simulation time of the last state update [s]
    scalar pocketTv_;                // vapour temperature (flux-weighted surface T) [K]
    // Local mode: one vapour column per radial bin (no radial vapour flow)
    scalarList pocketBinMass_;       // vapour mass per bin [kg]
    scalarList pocketBinDisp_;       // film displacement away from the glass per bin [m]
    scalarList pocketBinTv_;         // vapour temperature per bin [K]
    word pocketSavedMode_;           // pocket mode stored in the restart state ("" if none)
    scalar pocketSavedBinWidth_;     // bin width stored in the restart state (0 if none)

public:
    // Constructor
    advancedInterfaceCapturing
    (
        const fvMesh& mesh,
        volScalarField& alpha1,
        const twoPhaseMixtureThermo& mixture,
        const volScalarField& T,
        const volScalarField* Tl = nullptr
    );
    //- Delete default copy constructor and assignment operator
    advancedInterfaceCapturing(const advancedInterfaceCapturing&) = delete;
    void operator=(const advancedInterfaceCapturing&) = delete;
    // Member functions
    void correct();  // Main solver routine
    void calculateRecoilPressure();  // Calculate recoil pressure
    // Access functions
    const volScalarField& recoilPressure() const { return recoilPressure_; }
    dimensionedScalar meltingTemperature() const { return meltingTemp_; }
    dimensionedScalar vaporTemperature() const { return vaporTemp_; }
    // Write functions
    void write() const;
    void clearInternalFields()
    {
        // Reset internal fields to prevent memory issues
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        previousRecoilPressure_.setSize(mesh_.nCells());
        previousRecoilPressure_ = 0.0;
        havePreviousRecoil_ = false;
        rampProgress_ = 0.0;
        relaxTimeIndex_ = -1;
    }
};
} // End namespace Foam
#endif
