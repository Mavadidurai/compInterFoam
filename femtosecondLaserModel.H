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
    Foam::femtosecondLaserModel

Description
    Laser energy deposition model for Laser-Induced Forward Transfer (LIFT)
    simulations, providing temporal envelope handling, spatial weighting, and
    diagnostic reporting for femtosecond laser pulses.
\*---------------------------------------------------------------------------*/
#ifndef femtosecondLaserModel_H
#define femtosecondLaserModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "mathematicalConstants.H"
#include "tmp.H"
#include "typeInfo.H"
namespace Foam
{

class femtosecondLaserModel
{
private:
    // Mesh & config
    const fvMesh&       mesh_;
    const dictionary&   dict_;

    // Laser parameters
    dimensionedScalar   peakIntensity_;     // [W/m^2]
    dimensionedScalar   pulseWidth_;        // [s]
    dimensionedScalar   wavelength_;        // [m]
    dimensionedScalar   absorptionCoeff_;   // [1/m]
    dimensionedScalar   gasAbsorptionCoeff_; // [1/m]
    dimensionedScalar   spotSize_;          // [m] (diameter)
    dimensionedScalar   pulseEnergy_;       // [J]
    dimensionedScalar   maxVolumetricSource_; // [W/m^3]

    // Geometry / timing
    vector              direction_;         // unit vector
    point               focus_;
    point               initialFocus_;
    vector              scanVelocity_;      // [m/s]
    scalar              pulseFrequency_;    // [1/s]
    scalar              pulseDutyCycle_;    // [-]
    // Reflectivity specified at the film interface; converted to transmission
    // before spatial weighting so it attenuates the beam entering any cell.
    scalar              reflectivity_;      // [-]
    scalar              transmission_;      // [-] (>=0 overrides reflectivity)
    scalar              incidenceAngle_;    // [rad]
    bool                gaussianProfile_;
    bool                continuousLaser_;
    bool                applyInterfaceTransmissionEverywhere_;
    
    // Pulse energy consistency tolerances
    scalar              pulseEnergyToleranceRel_;
    scalar              pulseEnergyToleranceAbs_;
    
    scalar              laserStartTime_;    // [s]
    scalar              laserEndTime_;      // [s]
    scalar              pulseCenterTime_;   // [s] optional Gaussian centre
    // Donor film bounds (y)
    scalar              filmYMin_;
    scalar              filmYMax_;

    // Cached source backing the const source() API; update() resets these so
    // they are not treated as thread-safe caches
    mutable tmp<volScalarField> tSource_;
    mutable bool         sourceValid_;
    mutable scalar       cumulativeEnergy_;
    mutable scalar       cumulativeIncidentEnergy_;
    mutable scalar       cumulativeFilmEnergy_;
    mutable scalar       cumulativeGasEnergy_;
    mutable label        lastTimeIndex_;
    mutable scalar       pulseEnergyAccumulator_;
    mutable scalar       pulseExpectedAccumulator_;
    mutable scalar       pulseDepositableExpectedAccumulator_;
    // Integral of the source field as computed from the Gaussian/absorption
    // model, accumulated BEFORE the energy-conservation rescale in
    // calculateSource() forces it to match depositableExpectedThisStep.
    // Comparing this against the expected energy is the only check that can
    // actually catch a wrong Gaussian prefactor, FWHM conversion, or
    // absorption-depth bug -- the rescaled/deposited comparison cannot,
    // since the rescale makes it match the expectation by construction.
    mutable scalar       pulseRawEnergyAccumulator_;
    mutable scalar       currentPulseStartTime_;
    mutable label        pulseCounter_;
    mutable bool         trackingPulse_;
    mutable bool         pulseCompleted_;
    mutable bool         activeThisStep_;

    // Temporal envelope summary for the current time-step
    struct EnvelopeResult
    {
        scalar temporalIntegral = 0.0; // [s] integrated pulse envelope
        scalar temporalAverage  = 0.0; // [-] envelope averaged over dt
        scalar expectedEnergy   = 0.0; // [J] energy implied by envelope
        bool   active           = false; // true when the laser is on
    };

    // Accumulated spatial metrics during source assembly
    struct SpatialMetrics
    {
        label  cellsInBeam             = 0;   // cells touched by the beam
        label  cellsInFilm             = 0;   // beam cells within the film
        label  cellsInGas              = 0;   // beam cells within the gas
        label  limitedCells            = 0;   // cells capped at max source
        scalar maxSourceValue          = 0.0; // [W/m^3] max volumetric source
        scalar totalSourceIntegral     = 0.0; // [W] integrated volumetric source
        scalar totalFilmSourceIntegral = 0.0; // [W] power into the film
        scalar totalGasSourceIntegral  = 0.0; // [W] power into the gas
        scalar totalBeamVolume         = 0.0; // [m^3] summed beam volume
        scalar maxSourceCap            = 0.0; // [W/m^3] configured cap
        bool   limitSource             = false; // true when limiting is active
        scalar rawSourceIntegral       = 0.0; // [W] integrated source before
                                               // the energy-conservation rescale
    };

    // --- helpers ---
    scalar pulseNormalizationFactor() const;
    void   calculateSource() const;
    void   validateParameters() const;
    bool   checkPhysicalBounds() const;
    scalar calculateGaussianIntensity(const scalar R) const;
    bool   isInBeam(const point& p, const scalar axialHalfLength) const;
    bool   checkEnergyConservation() const;
    void   finalizePulseEnergyCheck(const char* context, const scalar currentTime) const;
    volScalarField& resetSourceField() const;  // ensure Q_laser exists and is zeroed
    scalar depositableEnergyFraction() const;
    scalar beamCoverageFraction() const;
    scalar effectiveTransmission(const scalar reflectivity) const;
    
    EnvelopeResult evaluateTemporalEnvelope
    (
        const scalar overlapStart,
        const scalar overlapEnd,
        const scalar dt
    ) const;
    scalar gaussianWindowIntegral
    (
        const scalar a,
        const scalar b,
        const scalar center,
        const scalar sigma
    ) const;
    SpatialMetrics applySpatialWeighting
    (
        volScalarField& source,
        const scalar temporalAverage
    ) const;
    void updateEnergyTracking
    (
        const volScalarField& source,
        const dimensionedScalar& dtDim,
        const EnvelopeResult& envelope,
        const SpatialMetrics& metrics,
        const scalar overlapStart,
        const scalar currentTime,
        const bool depositionActive,
        const label timeIndex
    ) const;
    void emitDiagnostics
    (
        const EnvelopeResult& envelope,
        const SpatialMetrics& metrics,
        const dimensionedScalar& dtDim,
        const dimensionedScalar& energyThisStep,
        const scalar incidentEnergyThisStep,
        const scalar filmEnergyThisStep,
        const scalar gasEnergyThisStep,
        const scalar currentTime,
        const label timeIndex,
        const bool depositionActive
    ) const;

    // disallow copy
    femtosecondLaserModel(const femtosecondLaserModel&) = delete;
    void operator=(const femtosecondLaserModel&) = delete;

public:
    TypeName("femtosecondLaserModel");

    // ctor/dtor
    femtosecondLaserModel(const fvMesh& mesh, const dictionary& dict);
    virtual ~femtosecondLaserModel() = default;

    // evolution
    void update();                   // scan + invalidate cache
    virtual void correct()           // force recompute source
    {
        sourceValid_ = false;
        calculateSource();
    }
    virtual void correct(const scalar) { correct(); }

    // access
    tmp<volScalarField> source() const;
    bool activeThisStep() const;
    const dimensionedScalar& peakIntensity() const { return peakIntensity_; }
    const dimensionedScalar& pulseWidth()    const { return pulseWidth_; }
    const dimensionedScalar& wavelength()    const { return wavelength_; }
    const dimensionedScalar& spotSize()      const { return spotSize_; }
    const dimensionedScalar& pulseEnergy()   const { return pulseEnergy_; }
    const dimensionedScalar& maxVolumetricSource() const {        return maxVolumetricSource_;    }
    const dimensionedScalar& gasAbsorptionCoeff() const
    {
        return gasAbsorptionCoeff_;
    }  
    const point&              focus()        const { return focus_; }
    const vector&             direction()    const { return direction_; }
    scalar laserStartTime() const { return laserStartTime_; }
    scalar laserEndTime()   const { return laserEndTime_; }
    scalar pulseCenterTime() const { return pulseCenterTime_; }    
    bool   isContinuous()   const { return continuousLaser_; }
    bool valid() const;
    void write() const;
};

} // End namespace Foam
#endif

