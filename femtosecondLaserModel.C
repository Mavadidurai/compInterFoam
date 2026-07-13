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

Description
    Implementation of the femtosecond laser source model used for
    Laser-Induced Forward Transfer (LIFT) simulations, including detailed
    diagnostics and optional temporal/spatial weighting controls.
\*---------------------------------------------------------------------------*/
#include "femtosecondLaserModel.H"
#include "fvc.H"
#include "fvm.H"
#include "mathematicalConstants.H"
#include "boundBox.H"
#include "treeBoundBox.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include <cmath>

extern Foam::Switch verbose;

namespace Foam
{
defineTypeNameAndDebug(femtosecondLaserModel, 0);

//------------------------------------------------------------------------------
// ctor
//------------------------------------------------------------------------------
femtosecondLaserModel::femtosecondLaserModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    peakIntensity_("peakIntensity", dimPower/dimArea,
        dict.getOrDefault<scalar>("peakIntensity", 0.0)),
    pulseWidth_("pulseWidth", dimTime, dict.get<scalar>("pulseWidth")),
    wavelength_("wavelength", dimLength, dict.get<scalar>("wavelength")),
    absorptionCoeff_("absorptionCoeff", dimless/dimLength,
    dict.getOrDefault<scalar>("absorptionCoeff", 5e6)),
    gasAbsorptionCoeff_(
        "gasAbsorptionCoeff",
        dimless/dimLength,
        dict.getOrDefault<scalar>("gasAbsorptionCoeff", 0.0)
    ),
    spotSize_("spotSize", dimLength, dict.get<scalar>("spotSize")),
    pulseEnergy_("pulseEnergy", dimEnergy, dict.get<scalar>("pulseEnergy")),
    maxVolumetricSource_
    (
        "maxVolumetricSource",
        dimPower/dimVolume,
        Foam::max
        (
            scalar(0),
            dict.getOrDefault<scalar>("maxVolumetricSource", 0.0)
        )
    ),
    direction_(dict.get<vector>("direction")),
    focus_(dict.get<point>("focus")),
    initialFocus_(dict.getOrDefault<point>("initialFocus", focus_)),
    scanVelocity_(dict.getOrDefault<vector>("scanVelocity", vector::zero)),
    pulseFrequency_(dict.getOrDefault<scalar>("pulseFrequency", 0.0)),
    pulseDutyCycle_(dict.getOrDefault<scalar>("pulseDutyCycle", 1.0)),
    reflectivity_(dict.getOrDefault<scalar>("reflectivity", 0.05)),
    transmission_(dict.getOrDefault<scalar>("transmission", -1.0)),
    incidenceAngle_(dict.getOrDefault<scalar>("incidenceAngle", 0.0)),
    gaussianProfile_(dict.getOrDefault<bool>("gaussianProfile", true)),
    continuousLaser_(dict.getOrDefault<bool>("continuousLaser", false)),
    applyInterfaceTransmissionEverywhere_
    (
        dict.getOrDefault<bool>("applyInterfaceTransmissionEverywhere", true)
    ),    
    pulseEnergyToleranceRel_
    (
        dict.getOrDefault<scalar>("pulseEnergyToleranceRel", 0.01)
    ),
    pulseEnergyToleranceAbs_
    (
        dict.getOrDefault<scalar>("pulseEnergyToleranceAbs", 1e-12)
    ),
    laserStartTime_(dict.getOrDefault<scalar>("laserStartTime", 0.0)),
    laserEndTime_(dict.get<scalar>("laserEndTime")),
    pulseCenterTime_(-GREAT),
    filmYMin_(0.0),
    filmYMax_(0.0),
    tSource_(),
    sourceValid_(false),
    cumulativeEnergy_(0.0),
    cumulativeIncidentEnergy_(0.0),
    cumulativeFilmEnergy_(0.0),
    cumulativeGasEnergy_(0.0),
    lastTimeIndex_(mesh.time().timeIndex()),
    pulseEnergyAccumulator_(0.0),
    pulseExpectedAccumulator_(0.0),
    pulseDepositableExpectedAccumulator_(0.0),
    pulseRawEnergyAccumulator_(0.0),
    currentPulseStartTime_(0.0),
    pulseCounter_(0),
    trackingPulse_(false),
    pulseCompleted_(false),
    activeThisStep_(false)
{
    const scalar originalReflectivity = reflectivity_;
    const scalar clampedReflectivity =
        Foam::min(Foam::max(originalReflectivity, scalar(0)), scalar(1));

    if (clampedReflectivity != originalReflectivity)
    {
        WarningInFunction
            << "reflectivity " << originalReflectivity
            << " outside [0, 1]; clamping to " << clampedReflectivity
            << endl;
    }

    reflectivity_ = clampedReflectivity;
    if (transmission_ >= 0)
    {
        const scalar originalTransmission = transmission_;
        const scalar clampedTransmission =
            Foam::min(Foam::max(originalTransmission, scalar(0)), scalar(1));

        if (clampedTransmission != originalTransmission)
        {
            WarningInFunction
                << "transmission " << originalTransmission
                << " outside [0, 1]; clamping to " << clampedTransmission
                << endl;
        }

        transmission_ = clampedTransmission;

        Info<< "Transmission override active: using transmission = "
            << transmission_ << "; reflectivity entry ignored" << endl;
    }

    if (pulseEnergyToleranceRel_ < 0)
    {
        WarningInFunction
            << "pulseEnergyToleranceRel " << pulseEnergyToleranceRel_
            << " is negative; clamping to 0" << endl;
        pulseEnergyToleranceRel_ = 0.0;
    }

    if (pulseEnergyToleranceAbs_ < 0)
    {
        WarningInFunction
            << "pulseEnergyToleranceAbs " << pulseEnergyToleranceAbs_
            << " is negative; clamping to 0" << endl;
        pulseEnergyToleranceAbs_ = 0.0;
    }

    if (dict.found("pulseCenterTime"))
    {
        pulseCenterTime_ = dict.get<scalar>("pulseCenterTime");
        const scalar clampedCenter =
            Foam::min(Foam::max(pulseCenterTime_, laserStartTime_), laserEndTime_);

        if (clampedCenter != pulseCenterTime_)
        {
            WarningInFunction
                << "pulseCenterTime " << pulseCenterTime_
                << " outside [" << laserStartTime_ << ", " << laserEndTime_
                << "]; clamping to " << clampedCenter
                << endl;
        }

        pulseCenterTime_ = clampedCenter;
    }
    const bool filmYMinProvided = dict.found("filmYMin");
    const bool filmYMaxProvided = dict.found("filmYMax");

    if (filmYMinProvided)
    {
        filmYMin_ = dict.get<scalar>("filmYMin");
    }

    if (filmYMaxProvided)
    {
        filmYMax_ = dict.get<scalar>("filmYMax");
    }

    if (filmYMinProvided != filmYMaxProvided)
    {
        FatalIOErrorInFunction(dict)
            << "Both 'filmYMin' and 'filmYMax' must be provided together"
            << " or omitted"
            << exit(FatalIOError);
    }

    if (!filmYMinProvided && !filmYMaxProvided)
    {
        scalar filmThickness = -1.0;
        word thicknessKey;

        if (dict.found("filmThickness"))
        {
            thicknessKey = "filmThickness";
            filmThickness = dict.get<scalar>(thicknessKey);
        }
        else if (dict.found("filmThicknessExpected"))
        {
            thicknessKey = "filmThicknessExpected";
            filmThickness = dict.get<scalar>(thicknessKey);
        }

        if (filmThickness <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "Specify positive film thickness via 'filmThickness'"
                << " (preferred) or 'filmThicknessExpected', or provide"
                << " explicit 'filmYMin'/'filmYMax' bounds"
                << exit(FatalIOError);
        }

        const scalar centerY = dict.lookupOrDefault<scalar>
        (
            "filmCenterY",
            initialFocus_.y()
        );
        const scalar halfThickness = 0.5*filmThickness;

        filmYMin_ = centerY - halfThickness;
        filmYMax_ = centerY + halfThickness;

        Info<< "Derived film bounds from center y=" << centerY
            << " m using " << thicknessKey << " = " << filmThickness
            << " m: [" << filmYMin_ << ", " << filmYMax_ << "] m" << endl;
    }
    else if (filmYMin_ >= filmYMax_)
    {
        FatalIOErrorInFunction(dict)
            << "Expected filmYMin < filmYMax; received "
            << filmYMin_ << " >= " << filmYMax_
            << exit(FatalIOError);
    }
    // normalize direction
    const scalar dirMag = mag(direction_);

    if (dirMag <= VSMALL)
    {
        FatalErrorInFunction
            << "Laser direction vector must be non-zero" << nl
            << "Supplied direction: " << direction_ << nl
            << abort(FatalError);
    }

    direction_ /= dirMag;
    if (pulseDutyCycle_ < 0.0 || pulseDutyCycle_ > 1.0)
    {
        const scalar originalDuty = pulseDutyCycle_;
        pulseDutyCycle_ = min(max(pulseDutyCycle_, 0.0), 1.0);

        WarningInFunction
            << "pulseDutyCycle " << originalDuty
            << " outside [0, 1]; clamping to " << pulseDutyCycle_ << endl;
    }
    // derive peak intensity if not provided
    const scalar normalization = max(VSMALL, pulseNormalizationFactor());
    const scalar derivedPeak = pulseEnergy_.value()/normalization;

    if (dict.found("peakIntensity"))
    {
        const scalar diff = mag(peakIntensity_.value() - derivedPeak);
        const scalar tol  = 0.05*derivedPeak;
        if (diff > tol)
        {
            const scalar adjustedEnergy = peakIntensity_.value()*normalization;
            const scalar originalEnergy = pulseEnergy_.value();

            WarningInFunction
                << "Supplied peakIntensity (" << peakIntensity_.value()
                << " W/m^2) is inconsistent with value derived from pulseEnergy "
                << "(" << derivedPeak << " W/m^2)." << nl
                << "Difference " << diff << " W/m^2 exceeds tolerance " << tol
                << " W/m^2; synchronizing configured pulse energy from "
                << originalEnergy << " J to " << adjustedEnergy
                << " J to match supplied intensity." << endl;

            pulseEnergy_ = dimensionedScalar
            (
                "pulseEnergy",
                dimEnergy,
                adjustedEnergy
            );
        }
    }
    else
    {
        peakIntensity_ = dimensionedScalar
        (
            "peakIntensity", dimPower/dimArea, derivedPeak
        );
    }

    if (!valid())
    {
        FatalErrorInFunction
            << "Invalid femtosecond laser configuration" << nl
            << "  Peak intensity: " << peakIntensity_.value() << nl
            << "  Pulse width: "    << pulseWidth_.value() << nl
            << "  Wavelength: "     << wavelength_.value() << nl
            << "  Spot size: "      << spotSize_.value() << nl
            << "  Pulse duty cycle: " << pulseDutyCycle_ << nl
            << abort(FatalError);
    }

    const bool master = Pstream::master();
    if (verbose && master)
    {
        Info<< "Femtosecond laser model initialized:" << nl
            << "  Mode: " << (continuousLaser_ ? "Continuous" : "Pulsed") << nl
            << "  Peak intensity: " << peakIntensity_.value() << " W/m^2" << nl
            << "  Pulse width: " << pulseWidth_.value() << " s" << nl
            << "  Wavelength: " << wavelength_.value() << " m" << nl
            << "  Spot size: " << spotSize_.value() << " m" << nl
            << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
            << "  Max volumetric source cap: " << maxVolumetricSource_.value()
            << " W/m^3" << nl
            << "  Focus: " << focus_ << nl
            << "  Direction: " << direction_ << nl
            << "  Active time: " << laserStartTime_ << " to "
            << laserEndTime_ << " s" << endl;
    }

    if (focus_.y() < filmYMin_ || focus_.y() > filmYMax_)
    {
        WarningInFunction
            << "Focus y=" << focus_.y()*1e6 << " µm outside film ["
            << filmYMin_*1e6 << ", " << filmYMax_*1e6
            << "] µm; coupling may be reduced." << endl;
    }
}

//------------------------------------------------------------------------------
// update (scan only, invalidate cache)
//------------------------------------------------------------------------------
void femtosecondLaserModel::update()
{
    if (mag(scanVelocity_) > VSMALL)
    {
        const scalar tnow = mesh_.time().value();
        focus_ = initialFocus_ + scanVelocity_ * tnow;
    }
    sourceValid_ = false;
}

//------------------------------------------------------------------------------
// validate parameters
//------------------------------------------------------------------------------
void femtosecondLaserModel::validateParameters() const
{
    if (peakIntensity_.dimensions() != dimPower/dimArea)
    {
        FatalErrorInFunction
            << "peakIntensity dimensions " << peakIntensity_.dimensions()
            << " do not match expected " << (dimPower/dimArea) << nl
            << abort(FatalError);
    }
    if (pulseWidth_.dimensions() != dimTime)
    {
        FatalErrorInFunction
            << "pulseWidth dimensions " << pulseWidth_.dimensions()
            << " do not match expected " << dimTime << nl
            << abort(FatalError);
    }
    if (wavelength_.dimensions() != dimLength)
    {
        FatalErrorInFunction
            << "wavelength dimensions " << wavelength_.dimensions()
            << " do not match expected " << dimLength << nl
            << abort(FatalError);
    }
    if (absorptionCoeff_.dimensions() != dimless/dimLength)
    {
        FatalErrorInFunction
            << "absorptionCoeff dimensions " << absorptionCoeff_.dimensions()
            << " do not match expected " << (dimless/dimLength) << nl
            << abort(FatalError);
    }
    if (gasAbsorptionCoeff_.dimensions() != dimless/dimLength)
    {
        FatalErrorInFunction
            << "gasAbsorptionCoeff dimensions "
            << gasAbsorptionCoeff_.dimensions()
            << " do not match expected " << (dimless/dimLength) << nl
            << abort(FatalError);
    }
    if (spotSize_.dimensions() != dimLength)
    {
        FatalErrorInFunction
            << "spotSize dimensions " << spotSize_.dimensions()
            << " do not match expected " << dimLength << nl
            << abort(FatalError);
    }
    if (pulseEnergy_.dimensions() != dimEnergy)
    {
        FatalErrorInFunction
            << "pulseEnergy dimensions " << pulseEnergy_.dimensions()
            << " do not match expected " << dimEnergy << nl
            << abort(FatalError);
    }
    if (maxVolumetricSource_.dimensions() != dimPower/dimVolume)
    {
        FatalErrorInFunction
            << "maxVolumetricSource dimensions "
            << maxVolumetricSource_.dimensions()
            << " do not match expected " << (dimPower/dimVolume) << nl
            << abort(FatalError);
    }

    if (peakIntensity_.value() <= 0)
    {
        FatalErrorInFunction
            << "peakIntensity (" << peakIntensity_.value()
            << ") must be positive" << nl
            << abort(FatalError);
    }
    if (pulseWidth_.value() <= 0)
    {
        FatalErrorInFunction
            << "pulseWidth (" << pulseWidth_.value()
            << ") must be positive" << nl
            << abort(FatalError);
    }
    if (wavelength_.value() <= 0)
    {
        FatalErrorInFunction
            << "wavelength (" << wavelength_.value()
            << ") must be positive" << nl
            << abort(FatalError);
    }
    if (spotSize_.value() <= 0)
    {
        FatalErrorInFunction
            << "spotSize (" << spotSize_.value()
            << ") must be positive" << nl
            << abort(FatalError);
    }
    if (pulseEnergy_.value() <= 0)
    {
        FatalErrorInFunction
            << "pulseEnergy (" << pulseEnergy_.value()
            << ") must be positive" << nl
            << abort(FatalError);
    }
    if (gasAbsorptionCoeff_.value() < 0)
    {
        FatalErrorInFunction
            << "gasAbsorptionCoeff (" << gasAbsorptionCoeff_.value()
            << ") must be non-negative" << nl
            << abort(FatalError);
    }
    if (maxVolumetricSource_.value() < 0)
    {
        FatalErrorInFunction
            << "maxVolumetricSource (" << maxVolumetricSource_.value()
            << ") must be non-negative" << nl
            << abort(FatalError);
    }
    if (laserStartTime_ < 0)
    {
        FatalErrorInFunction
            << "laserStartTime (" << laserStartTime_
            << ") must be non-negative" << nl
            << abort(FatalError);
    }

    if (laserEndTime_ <= laserStartTime_)
    {
        FatalErrorInFunction
            << "laserEndTime (" << laserEndTime_
            << ") must be greater than laserStartTime ("
            << laserStartTime_ << ")" << nl
            << abort(FatalError);
    }
    if (pulseDutyCycle_ < 0.0 || pulseDutyCycle_ > 1.0)
    {
        FatalErrorInFunction
            << "pulseDutyCycle (" << pulseDutyCycle_
            << ") must be within [0, 1]" << nl
            << abort(FatalError);
    }

    // film thickness sanity
    const scalar filmThickness = filmYMax_ - filmYMin_;
    scalar expected = filmThickness;

    if (dict_.found("filmThickness"))
    {
        expected = dict_.get<scalar>("filmThickness");
    }
    else if (dict_.found("filmThicknessExpected"))
    {
        expected = dict_.get<scalar>("filmThicknessExpected");
    }

    const scalar tol =
        dict_.lookupOrDefault<scalar>
        (
            "filmThicknessTolerance",
            0.1*Foam::max(expected, scalar(SMALL))
        );
    if (filmThickness <= 0)
    {
        WarningInFunction
            << "Non-positive film thickness: " << filmThickness
            << " (check filmYMin, filmYMax)" << endl;
    }
    if (expected <= SMALL)
    {
        WarningInFunction
            << "Configured film thickness expectation " << expected
            << " m is non-positive" << endl;
    }
    else if (mag(filmThickness - expected) > tol)
    {
        WarningInFunction
            << "Donor film thickness (" << filmThickness
            << ") deviates from expected " << expected
            << " by more than " << tol << " m" << endl;
    }
}

//------------------------------------------------------------------------------
// physical bounds
//------------------------------------------------------------------------------
bool femtosecondLaserModel::checkPhysicalBounds() const
{
    bool ok = true;

    if (!continuousLaser_
     && (pulseWidth_.value() < 1e-16 || pulseWidth_.value() > 1e-11))
    {
        WarningInFunction << "Pulse width out of 10fs–10ps range: "
                          << pulseWidth_.value() << " s" << endl;
        //ok = false;
    }
    if (wavelength_.value() < 1e-7 || wavelength_.value() > 2e-6)
    {
        WarningInFunction << "Wavelength out of 100nm–2µm range: "
                          << wavelength_.value() << " m" << endl;
        //ok = false;
    }
    // Typical LIFT experiments routinely exceed 6e16 W/m^2, so only warn
    // when intensities move well beyond the validated range for this model.
    if (peakIntensity_.value() > 1e17)
    {
        WarningInFunction << "Very high peak intensity: "
                          << peakIntensity_.value() << " W/m^2" << endl;
        //ok = false;
    }
    return ok;
}
//------------------------------------------------------------------------------
scalar femtosecondLaserModel::pulseNormalizationFactor() const
{
    const scalar beamRadius = spotSize_.value()/2.0;
    scalar spatialFactor = constant::mathematical::pi*sqr(beamRadius);

    if (gaussianProfile_)
    {
        spatialFactor *= 0.5;
    }

    scalar temporalFactor = pulseWidth_.value();

    if (!continuousLaser_)
    {
        const scalar sqrtPi = sqrt(constant::mathematical::pi);
        const scalar sqrtLn2 = sqrt(log(2.0));
        temporalFactor *= sqrtPi/(2.0*sqrtLn2);
    }

    return spatialFactor*max(VSMALL, temporalFactor);
}
//------------------------------------------------------------------------------
scalar femtosecondLaserModel::calculateGaussianIntensity(const scalar R) const
{
    if (gaussianProfile_)
    {
        const scalar w = spotSize_.value()/2.0; // radius
        return exp(-2.0*sqr(R)/sqr(w));
    }
    else
    {
        return (R <= spotSize_.value()/2.0) ? 1.0 : 0.0;
    }
}

//------------------------------------------------------------------------------
bool femtosecondLaserModel::isInBeam
(
    const point& p,
    const scalar axialHalfLength
) const
{
    vector r = p - focus_;
    scalar z = (r & direction_);
    r -= z*direction_;
    const scalar R = mag(r);

    const scalar beamRadius = spotSize_.value()/2.0;
    const scalar maxRadius  = 3.0*beamRadius;         // ~3-sigma
    const scalar maxZ       = max(axialHalfLength, 0.0);

    return (R <= maxRadius) && (z >= -maxZ) && (z <= maxZ);
}

//------------------------------------------------------------------------------
bool femtosecondLaserModel::checkEnergyConservation() const
{
    if (!tSource_.valid() || continuousLaser_) return true;

    const dimensionedScalar dt = mesh_.time().deltaT();
    const dimensionedScalar totalEnergy = fvc::domainIntegrate(tSource_()*dt);

    const dimensionedScalar expected
    (
        "expectedEnergy", dimEnergy,
        peakIntensity_.value()*pulseNormalizationFactor()
    );
    const scalar maxAllowed = 10.0*expected.value(); // very generous
    return totalEnergy.value() <= maxAllowed;
}
//------------------------------------------------------------------------------
void femtosecondLaserModel::finalizePulseEnergyCheck
(
    const char* context,
    const scalar currentTime
) const
{
    if (continuousLaser_ || !trackingPulse_) return;

    if (pulseExpectedAccumulator_ <= VSMALL)
    {
        trackingPulse_ = false;
        pulseEnergyAccumulator_ = 0.0;
        pulseExpectedAccumulator_ = 0.0;
        pulseDepositableExpectedAccumulator_ = 0.0;
        pulseRawEnergyAccumulator_ = 0.0;
        currentPulseStartTime_ = 0.0;
        return;
    }

    const scalar incidentExpected = pulseExpectedAccumulator_;
    const scalar configuredIncident = pulseEnergy_.value();
    const scalar depositableExpected =
        pulseDepositableExpectedAccumulator_;
    scalar depositableRatio = 0.0;

    if (incidentExpected > VSMALL)
    {
        depositableRatio = depositableExpected/incidentExpected;
    }
    else
    {
        depositableRatio =
            depositableEnergyFraction()*beamCoverageFraction();
    }

    depositableRatio = min(max(depositableRatio, scalar(0)), scalar(1));
    const scalar configuredDepositable = configuredIncident*depositableRatio;

    // Do not artificially inflate the tolerance reference for sub-nJ pulses;
    // scale it with the actual expected/depositable energy so the audit
    // continues to react when very small pulses are under-delivered.
    const scalar reference =
        max(max(depositableExpected, configuredDepositable), scalar(VSMALL));
    const scalar tolerance =
        max(pulseEnergyToleranceAbs_, pulseEnergyToleranceRel_*reference);

    const scalar diffDepositableExpected =
        mag(pulseEnergyAccumulator_ - depositableExpected);
    const scalar diffConfigured =
        mag(pulseEnergyAccumulator_ - configuredDepositable);
    const scalar diffIncident =
        mag(pulseEnergyAccumulator_ - incidentExpected);

    const scalar configuredGap =
        mag(depositableExpected - configuredDepositable);
    const scalar incidentGap =
        mag(incidentExpected - depositableExpected);

    const scalar configuredResidual =
        max(diffConfigured - configuredGap, scalar(0));
    const scalar incidentResidual =
        max(diffIncident - incidentGap, scalar(0));
    const scalar maxDeviation =
        max(diffDepositableExpected, diffConfigured);

    // Independent check: pulseRawEnergyAccumulator_ is the integral of the
    // source field as computed by the Gaussian/absorption model BEFORE the
    // energy-conservation rescale above forces its total to match
    // depositableExpected. Unlike diffDepositableExpected (which compares
    // the rescaled/deposited energy and is tautologically close to zero by
    // construction), this can actually catch a wrong Gaussian prefactor,
    // FWHM conversion, or absorption-depth bug in the spatial model.
    const scalar diffRawModel =
        mag(pulseRawEnergyAccumulator_ - depositableExpected);
    const bool warnRawModel = diffRawModel > tolerance;

    ++pulseCounter_;

    const bool master = Pstream::master();
    if (verbose && master)
    {
        Info<< "LASER PULSE ENERGY CHECK:" << nl
            << "  Pulse index:      " << pulseCounter_ << nl
            << "  Context:          " << context << nl
            << "  Start time:       " << currentPulseStartTime_ << " s" << nl
            << "  End time:         " << currentTime << " s" << nl
            << "  Deposited energy:           " << pulseEnergyAccumulator_
            << " J" << nl
            << "  Raw model energy (pre-rescale): " << pulseRawEnergyAccumulator_
            << " J" << nl
            << "  Incident expected (integrated): " << incidentExpected << " J" << nl
            << "  Depositable expected (integrated): " << depositableExpected
            << " J" << nl
            << "  Depositable/incident ratio: " << depositableRatio << nl
            << "  Configured incident:        " << configuredIncident << " J" << nl
            << "  Configured depositable:     " << configuredDepositable << " J"
            << nl
            << "  Diff vs depositable expected: " << diffDepositableExpected
            << " J" << nl
            << "  Diff vs configured depositable:   " << diffConfigured
            << " J" << nl
            << "  Diff vs incident expected:   " << diffIncident << " J" << nl
            << "  Configured depositable gap:  " << configuredGap << " J" << nl
            << "  Configured depositable residual:  " << configuredResidual
            << " J" << nl
            << "  Incident gap:          " << incidentGap << " J" << nl
            << "  Incident residual:     " << incidentResidual << " J" << nl
            << "  Max deviation:        " << maxDeviation << " J" << nl
            << "  Tolerance(abs):       " << pulseEnergyToleranceAbs_ << " J" << nl
            << "  Tolerance(rel*ref):   "
            << pulseEnergyToleranceRel_*reference << " J" << nl
            << "  Applied tolerance:    " << tolerance << " J" << endl;
    }

    const bool warnDepositable = diffDepositableExpected > tolerance;
    const bool warnConfigured = configuredResidual > tolerance;
    const bool warnIncident = incidentResidual > tolerance;

    if (warnRawModel)
    {
        WarningInFunction
            << "Laser spatial source model mismatch after pulse "
            << pulseCounter_ << ": the Gaussian/absorption model integrated "
            << "to " << pulseRawEnergyAccumulator_
            << " J before energy-conservation rescaling, vs depositable "
            << "energy expected from the temporal envelope "
            << depositableExpected << " J. Residual mismatch "
            << diffRawModel << " J exceeds tolerance " << tolerance
            << " J -- check the Gaussian prefactor, FWHM conversion, and"
            << " absorption depth; the rescale step masks this from the"
            << " deposited-energy check above." << endl;
    }

    if (warnDepositable)
    {
        WarningInFunction
            << "Laser pulse energy mismatch after pulse " << pulseCounter_
            << ": deposited " << pulseEnergyAccumulator_
            << " J vs depositable energy expected from temporal envelope "
            << depositableExpected << " J. Residual mismatch "
            << diffDepositableExpected << " J exceeds tolerance "
            << tolerance << " J" << endl;
    }

    if (warnConfigured)
    {
        WarningInFunction
            << "Laser pulse energy mismatch after pulse " << pulseCounter_
            << ": deposited " << pulseEnergyAccumulator_
            << " J vs configured depositable energy " << configuredDepositable
            << " J (reflective/transmissive losses excluded). Residual mismatch "
            << configuredResidual << " J exceeds tolerance "
            << tolerance << " J" << endl;
    }

    if (warnIncident)
    {
        WarningInFunction
            << "Laser pulse energy mismatch after pulse " << pulseCounter_
            << ": deposited " << pulseEnergyAccumulator_ << " J vs incident"
            << " energy expectation " << incidentExpected
            << " J after accounting for interface/coverage losses (residual "
            << incidentResidual << " J exceeds tolerance "
            << tolerance << " J)." << endl;
    }

    if (!warnDepositable && !warnConfigured && !warnIncident && verbose && master)
    {
        if (diffDepositableExpected > VSMALL)
        {
            Info<< "  Note: deposited energy differs from depositable expectation"
                << " by " << diffDepositableExpected
                << " J (within tolerance)." << endl;
        }

        if (diffConfigured > VSMALL)
        {
            Info<< "  Note: deposited energy differs from configured depositable"
                << " energy by " << diffConfigured
                << " J (within tolerance)." << endl;
        }

        if (diffIncident > VSMALL)
        {
            Info<< "  Note: deposited energy differs from incident expectation"
                << " by " << diffIncident << " J (within tolerance)." << endl;
        }

        if (configuredGap > VSMALL && configuredGap <= tolerance)
        {
            Info<< "  Note: depositable expectation differs from configured"
                << " depositable energy by " << configuredGap
                << " J (within tolerance)." << endl;
        }
    }

    trackingPulse_ = false;
    pulseEnergyAccumulator_ = 0.0;
    pulseExpectedAccumulator_ = 0.0;
    pulseDepositableExpectedAccumulator_ = 0.0;
    pulseRawEnergyAccumulator_ = 0.0;
    currentPulseStartTime_ = 0.0;
}

//------------------------------------------------------------------------------
volScalarField& femtosecondLaserModel::resetSourceField() const
{
    if (!tSource_.valid())
    {
        tSource_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "Q_laser",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimPower/dimVolume, 0.0)
            )
        );
    }

    volScalarField& source = tSource_.ref();
    source = dimensionedScalar("zero", source.dimensions(), 0.0);
    return source;
}

//------------------------------------------------------------------------------
// Integrate a Gaussian pulse over [a, b] for energy bookkeeping
//------------------------------------------------------------------------------
scalar femtosecondLaserModel::gaussianWindowIntegral
(
    const scalar a,
    const scalar b,
    const scalar center,
    const scalar sigma
) const
{
    if (sigma <= 0)
    {
        return 0.0;
    }

    const scalar invSqrt2Sigma = 1.0/(sqrt(2.0)*sigma);
    const scalar prefactor = sigma*sqrt(constant::mathematical::pi/2.0);

    return prefactor
        * (std::erf((b - center)*invSqrt2Sigma)
         - std::erf((a - center)*invSqrt2Sigma));
}
//------------------------------------------------------------------------------
scalar femtosecondLaserModel::effectiveTransmission(const scalar reflectivity) const
{
    scalar transmissionFactor = 1.0 - reflectivity;

    if (transmission_ >= 0)
    {
        transmissionFactor = transmission_;
    }
    else if (dict_.found("reflectivity"))
    {
        transmissionFactor = 1.0 - reflectivity;
    }
    else
    {
        // Fresnel reflection at the air-metal interface (normal incidence)
        // Empirical optical constants for Ti at 343 nm
        const scalar n_metal = 2.0;
        const scalar k_metal = 2.8;
        const scalar n_air = 1.0;

        const scalar num_real = n_air - n_metal;
        const scalar num_imag = -k_metal;
        const scalar den_real = n_air + n_metal;
        const scalar den_imag = k_metal;

        const scalar R_fresnel =
            (sqr(num_real) + sqr(num_imag)) /
            (sqr(den_real) + sqr(den_imag));

        const scalar T_interface = 1.0 - R_fresnel;

        const scalar filmThickness = Foam::max(filmYMax_ - filmYMin_, scalar(0));
        const scalar alpha_abs = absorptionCoeff_.value();
        const scalar opticalDepth = alpha_abs*filmThickness;
        const scalar R_back = R_fresnel;
        const scalar denom = 1.0 - R_back*R_back*std::exp(-2.0*opticalDepth);
        const scalar multipleReflFactor =
            (mag(denom) > VSMALL) ? (1.0/denom) : 1.0;

        transmissionFactor = T_interface * multipleReflFactor;
    }

    return min(max(transmissionFactor, scalar(0)), scalar(1));
}

scalar femtosecondLaserModel::depositableEnergyFraction() const
{
    const scalar transmissionFactor = effectiveTransmission(reflectivity_);

    const scalar filmThickness = max(filmYMax_ - filmYMin_, scalar(0));
    const scalar dirYMag = mag(direction_.y());
    scalar filmPathLength = 0.0;

    if (filmThickness > VSMALL)
    {
        if (dirYMag > VSMALL)
        {
            filmPathLength = filmThickness/dirYMag;
        }
        else
        {
            filmPathLength = filmThickness;
        }
    }

    auto absorptionFraction = [](const scalar opticalDepth)
    {
        if (opticalDepth <= VSMALL)
        {
            return scalar(0.0);
        }

        if (opticalDepth < 1e-6)
        {
            return opticalDepth - 0.5*sqr(opticalDepth);
        }

        return scalar(1.0 - Foam::exp(-opticalDepth));
    };

    scalar filmFraction = absorptionFraction
    (
        absorptionCoeff_.value()*filmPathLength
    );
    filmFraction = min(max(filmFraction, scalar(0)), scalar(1));

    scalar gasFraction = 0.0;
    const scalar gasCoeff = gasAbsorptionCoeff_.value();

    if (gasCoeff > VSMALL)
    {
        scalar gasDistance = 0.0;
        const boundBox& bounds = mesh_.bounds();
        const scalar dirY = direction_.y();

        if (dirY < -VSMALL)
        {
            gasDistance = max(bounds.max().y() - filmYMax_, scalar(0));
        }
        else if (dirY > VSMALL)
        {
            gasDistance = max(filmYMin_ - bounds.min().y(), scalar(0));
        }

        if (gasDistance > VSMALL && dirYMag > VSMALL)
        {
            gasDistance /= dirYMag;
        }
        else
        {
            gasDistance = 0.0;
        }

        gasFraction = absorptionFraction(gasCoeff*gasDistance);
        gasFraction = min(max(gasFraction, scalar(0)), scalar(1));
    }

    scalar netFraction = 0.0;

    if (transmissionFactor > VSMALL)
    {
        const scalar gasDepositable = transmissionFactor*gasFraction;
        const scalar filmDepositable =
            transmissionFactor*(1.0 - gasFraction)*filmFraction;

        netFraction = gasDepositable + filmDepositable;
    }

    return min(max(netFraction, scalar(0)), scalar(1));
}
//------------------------------------------------------------------------------
scalar femtosecondLaserModel::beamCoverageFraction() const
{
    if (!gaussianProfile_)
    {
        return 1.0;
    }

    const scalar beamRadius = spotSize_.value()/2.0;
    if (beamRadius <= VSMALL)
    {
        return 0.0;
    }

    const boundBox bounds = mesh_.bounds();
    if (!bounds.valid())
    {
        return 1.0;
    }

    vector refAxis(1, 0, 0);
    if (mag(direction_ ^ refAxis) <= VSMALL)
    {
        refAxis = vector(0, 1, 0);
    }
    if (mag(direction_ ^ refAxis) <= VSMALL)
    {
        refAxis = vector(0, 0, 1);
    }

    vector e1 = direction_ ^ refAxis;
    scalar e1Mag = mag(e1);
    if (e1Mag <= VSMALL)
    {
        return 1.0;
    }
    e1 /= e1Mag;

    vector e2 = direction_ ^ e1;
    const scalar e2Mag = mag(e2);
    if (e2Mag <= VSMALL)
    {
        return 1.0;
    }
    e2 /= e2Mag;

    const point minCorner = bounds.min();
    const point maxCorner = bounds.max();

    scalar minU = GREAT;
    scalar maxU = -GREAT;
    scalar minV = GREAT;
    scalar maxV = -GREAT;

    for (label i = 0; i < 8; ++i)
    {
        const point corner
        (
            (i & 1) ? maxCorner.x() : minCorner.x(),
            (i & 2) ? maxCorner.y() : minCorner.y(),
            (i & 4) ? maxCorner.z() : minCorner.z()
        );

        const vector rel = corner - focus_;
        const scalar u = rel & e1;
        const scalar v = rel & e2;

        minU = min(minU, u);
        maxU = max(maxU, u);
        minV = min(minV, v);
        maxV = max(maxV, v);
    }

    if (maxU <= minU || maxV <= minV)
    {
        return 0.0;
    }

    const scalar sqrt2OverW = sqrt(2.0)/beamRadius;
    const scalar fracU = 0.5
        * (std::erf(sqrt2OverW*maxU) - std::erf(sqrt2OverW*minU));
    const scalar fracV = 0.5
        * (std::erf(sqrt2OverW*maxV) - std::erf(sqrt2OverW*minV));

    scalar coverage = fracU*fracV;
    coverage = min(max(coverage, scalar(0)), scalar(1));

    return coverage; // ignores tree-box clipping; focus near bounds overestimates
}
//------------------------------------------------------------------------------
femtosecondLaserModel::EnvelopeResult
femtosecondLaserModel::evaluateTemporalEnvelope
(
    const scalar overlapStart,
    const scalar overlapEnd,
    const scalar dt
) const
{
    EnvelopeResult result;
    if (continuousLaser_)
    {
        result.temporalIntegral = overlapEnd - overlapStart;

        if (result.temporalIntegral > VSMALL)
        {
            const scalar beamRadius = spotSize_.value()/2.0;
            if (beamRadius > VSMALL)
            {
                scalar effectiveArea = constant::mathematical::pi*sqr(beamRadius);

                if (gaussianProfile_)
                {
                    effectiveArea *= 0.5;
                }

                const scalar incidentPower =
                    peakIntensity_.value()*effectiveArea;

                result.expectedEnergy =
                    incidentPower
                  * result.temporalIntegral;
            }
        }
    }
    else if (pulseFrequency_ > SMALL)
    {
        const scalar period   = 1.0/pulseFrequency_;
        const scalar onTime   = max(SMALL, pulseDutyCycle_*period);
        const scalar sigma    = pulseWidth_.value()
            /(2.0*sqrt(2.0*log(2.0)));
        const scalar gaussianNormalization =
            sigma*sqrt(constant::mathematical::twoPi);
        const scalar pulseWidthInv =
            sqrt(4.0*log(2.0))/max(pulseWidth_.value(), VSMALL);
        const scalar localStart = overlapStart - laserStartTime_;
        const scalar localEnd   = overlapEnd   - laserStartTime_;

        label periodIndex = 0;
        if (localStart > SMALL)
        {
            periodIndex = static_cast<label>(std::floor(localStart/period));
        }

        for
        (
            scalar windowStart = periodIndex*period;
            windowStart < localEnd + period;
            windowStart += period
        )
        {
            const scalar windowEnd = windowStart + onTime;

            if (windowEnd <= localStart)
            {
                continue;
            }
            if (windowStart >= localEnd)
            {
                break;
            }

            const scalar clipStart = max(windowStart, localStart);
            const scalar clipEnd   = min(windowEnd, localEnd);

            if (clipEnd <= clipStart)
            {
                continue;
            }

            const scalar center = windowStart + 0.5*onTime;
            const scalar argStart =
                (clipStart - center)*pulseWidthInv;
            const scalar argEnd =
                (clipEnd - center)*pulseWidthInv;
            const scalar fraction = 0.5
                * (std::erf(argEnd) - std::erf(argStart));

            if (mag(fraction) > VSMALL)
            {
                result.temporalIntegral += fraction*gaussianNormalization;
                result.expectedEnergy += pulseEnergy_.value()*fraction;
            }
        }
    }
    else
    {
        // Single Gaussian pulse.  Users often keep laserEndTime well beyond the
        // actual pulse duration so that other physics (e.g. phase change) remain
        // enabled.  The previous implementation centred the Gaussian at the
        // midpoint between laserStartTime and laserEndTime, which can shift the
        // peak several picoseconds after the intended trigger and effectively
        // delay the heat source.  Constrain the peak to remain close to the        
        // configured start time by capping how far it may drift from
        // laserStartTime when no explicit pulseCenterTime is supplied.

        const scalar sigma =
            pulseWidth_.value()/(2.0*sqrt(2.0*log(2.0)));
            // std dev from FWHM
        scalar center = pulseCenterTime_;
        const bool customCenter = (pulseCenterTime_ > -0.5*GREAT);

        if (!customCenter)
        {
            const scalar windowWidth = max(laserEndTime_ - laserStartTime_, SMALL);

            if (sigma <= VSMALL)
            {
                center = laserStartTime_ + 0.5*windowWidth;
            }
            else
            {
                const scalar halfWindow = 0.5*windowWidth;
                const scalar maxCenterOffset = 3.0*sigma;
                const scalar centerOffset = min(halfWindow, maxCenterOffset);

                center = laserStartTime_ + centerOffset;
            }
        }

        result.temporalIntegral =
            gaussianWindowIntegral(overlapStart, overlapEnd, center, sigma);

        if (result.temporalIntegral > VSMALL)
        {
            const scalar windowIntegral =
                gaussianWindowIntegral
                (
                    laserStartTime_,
                    laserEndTime_,
                    center,
                    sigma
                );
            result.expectedEnergy =
                pulseEnergy_.value()
              * result.temporalIntegral/max(windowIntegral, VSMALL);
        }
    }

    result.active = result.temporalIntegral > VSMALL;

    if (result.active)
    {
        result.temporalAverage = result.temporalIntegral/max(dt, VSMALL);
        result.temporalAverage =
            min(scalar(1.0), max(result.temporalAverage, scalar(0.0)));
    }

    return result;
}

//------------------------------------------------------------------------------
femtosecondLaserModel::SpatialMetrics
femtosecondLaserModel::applySpatialWeighting
(
    volScalarField& source,
    const scalar temporalAverage
) const
{
    SpatialMetrics metrics;
    metrics.maxSourceCap = maxVolumetricSource_.value();
    metrics.limitSource = metrics.maxSourceCap > SMALL;

    const scalar interfaceTransmission = effectiveTransmission(reflectivity_);
    const bool transmissionOverride = transmission_ >= 0;
    const bool reflectivityConfigured = dict_.found("reflectivity");

   if (verbose && Pstream::master())
    {
        Info<< "Interface transmission factor applied = "
            << interfaceTransmission;

        if (transmissionOverride)
        {
            Info<< " (user transmission override)";
        }
        else if (reflectivityConfigured)
        {
            Info<< " (from reflectivity = " << reflectivity_ << ")";
        }
        else
        {
            Info<< " (default Fresnel interface model)";
        }

        Info<< endl;
    }

    const pointField& cellCentres = mesh_.C();
    const scalarField& cellVolumes = mesh_.V();
    const pointField& meshPoints = mesh_.points();
    const volScalarField* metalFractionPtr = nullptr;

    if (mesh_.foundObject<volScalarField>("alpha.metal"))
    {
        metalFractionPtr = &mesh_.lookupObject<volScalarField>("alpha.metal");
    }
    else if (mesh_.foundObject<volScalarField>("alpha1"))
    {
        metalFractionPtr = &mesh_.lookupObject<volScalarField>("alpha1");
    }

    const scalarField* metalFractions = nullptr;

    if (metalFractionPtr)
    {
        metalFractions = &metalFractionPtr->internalField();
    }

    const scalar metalFractionTol = 1e-6;

    // filmYMin_/filmYMax_ are set once from the initial geometry and never
    // updated as the film melts, spreads, or ejects material. Rebuild the
    // effective band from the live alpha field each call so absorption
    // depth and ray-clipping below track the metal that is actually there
    // instead of silently zeroing deposition once it moves outside the
    // original band.
    scalar liveFilmYMin = filmYMin_;
    scalar liveFilmYMax = filmYMax_;

    if (metalFractions)
    {
        scalar localMin = GREAT;
        scalar localMax = -GREAT;
        bool foundMetal = false;

        forAll(*metalFractions, cellI)
        {
            if ((*metalFractions)[cellI] > metalFractionTol)
            {
                const scalar y = cellCentres[cellI].y();
                localMin = Foam::min(localMin, y);
                localMax = Foam::max(localMax, y);
                foundMetal = true;
            }
        }

        reduce(foundMetal, orOp<bool>());

        if (foundMetal)
        {
            reduce(localMin, minOp<scalar>());
            reduce(localMax, maxOp<scalar>());
            liveFilmYMin = localMin;
            liveFilmYMax = localMax;
        }
    }

    const scalar filmThickness = Foam::max(liveFilmYMax - liveFilmYMin, scalar(0));
    const scalar beamRadius = spotSize_.value()/2.0;
    const scalar radialHalfWidth = 3.0*beamRadius; // ~3-sigma laterally

    const scalar filmHalfThickness = 0.5*filmThickness;

    scalar axialHalfLength = Foam::max(filmHalfThickness, beamRadius);

    if (absorptionCoeff_.value() > VSMALL)
    {
        const scalar absorptionDepth = 1.0/absorptionCoeff_.value();
        axialHalfLength = Foam::max(axialHalfLength, 3.0*absorptionDepth);
    }

    axialHalfLength = Foam::max(axialHalfLength, SMALL);

    const vector axialContribution
    (
        axialHalfLength*mag(direction_.x()),
        axialHalfLength*mag(direction_.y()),
        axialHalfLength*mag(direction_.z())
    );

    const vector halfWidths =
        vector(radialHalfWidth, radialHalfWidth, radialHalfWidth)
      + axialContribution;

    const treeBoundBox searchBox(focus_ - halfWidths, focus_ + halfWidths);

    if (verbose && Pstream::master())
    {

        Info<< "==== LASER BEAM DIAGNOSTIC ====" << nl
            << "Focus position: " << focus_ << nl
            << "Beam direction: " << direction_ << nl
            << "Spot size (diameter): " << spotSize_.value() << " m" << nl
            << "Beam radius: " << beamRadius << " m" << nl
            << "Radial half-width (3σ): " << radialHalfWidth << " m" << nl
            << "Axial half-length: " << axialHalfLength << " m" << nl
            << "Search box min: " << searchBox.min() << nl
            << "Search box max: " << searchBox.max() << nl
            << "Mesh bounds min: " << mesh_.bounds().min() << nl
            << "Mesh bounds max: " << mesh_.bounds().max() << nl
            << "Total mesh cells: " << mesh_.nCells() << nl
            << "==============================" << endl;
    }

    const vector directionUnit(direction_);

    point entryPoint = focus_;

    const boundBox domainBox(mesh_.bounds());
    const point domainMin = domainBox.min();
    const point domainMax = domainBox.max();

    scalar tMin = -GREAT;
    scalar tMax = GREAT;

    for (Foam::direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const scalar dirComp = directionUnit[cmpt];
        const scalar focusComp = focus_[cmpt];

        if (mag(dirComp) < VSMALL)
        {
            if
            (
                focusComp < domainMin[cmpt] - SMALL
             || focusComp > domainMax[cmpt] + SMALL
            )
            {
                tMin = GREAT;
                tMax = -GREAT;
                break;
            }

            continue;
        }

        const scalar t1 = (domainMin[cmpt] - focusComp)/dirComp;
        const scalar t2 = (domainMax[cmpt] - focusComp)/dirComp;

        const scalar axisMin = Foam::min(t1, t2);
        const scalar axisMax = Foam::max(t1, t2);

        tMin = Foam::max(tMin, axisMin);
        tMax = Foam::min(tMax, axisMax);
    }

    if (tMin <= tMax && tMax >= 0)
    {
        const scalar entryParam =
            (tMin > 0) ? tMin : Foam::min(tMin, scalar(0));
        entryPoint = focus_ + entryParam*directionUnit;
    }
    scalar filmIntervalStart = 0.0;
    scalar filmIntervalEnd = 0.0;
    bool haveFilmInterval = false;

    const scalar dirY = directionUnit.y();

    if (mag(dirY) > VSMALL)
    {
        const scalar sToMin = (liveFilmYMin - entryPoint.y())/dirY;
        const scalar sToMax = (liveFilmYMax - entryPoint.y())/dirY;
        const scalar intervalMin = Foam::min(sToMin, sToMax);
        const scalar intervalMax = Foam::max(sToMin, sToMax);
        const scalar positiveStart = Foam::max(intervalMin, scalar(0));

        if (intervalMax > positiveStart)
        {
            filmIntervalStart = positiveStart;
            filmIntervalEnd = intervalMax;
            haveFilmInterval = true;
        }
    }
    else if
    (
        entryPoint.y() >= liveFilmYMin
     && entryPoint.y() <= liveFilmYMax
    )
    {
        filmIntervalStart = 0.0;
        filmIntervalEnd = GREAT;
        haveFilmInterval = true;
    }
    scalar filmEntryOffset = 0.0;

    if (haveFilmInterval)
    {
        const scalar dirY = directionUnit.y();

        if (mag(dirY) > VSMALL)
        {
            const scalar targetY = (dirY >= 0) ? liveFilmYMin : liveFilmYMax;
            filmEntryOffset = (targetY - entryPoint.y())/dirY;
            filmEntryOffset = Foam::max(filmEntryOffset, scalar(0));
        }
    }
    label cellsChecked = 0;
    label cellsInSearchBox = 0;
    
    forAll(cellCentres, cellI)
    {
        const point& c = cellCentres[cellI];

        if (!searchBox.contains(c))
        {
            continue;
        }

        ++cellsInSearchBox;
        ++cellsChecked;
        
        bool inFilm = (c.y() >= liveFilmYMin && c.y() <= liveFilmYMax);

        if (metalFractions)
        {
            const scalar alphaVal = Foam::min
            (
                Foam::max((*metalFractions)[cellI], scalar(0)),
                scalar(1)
            );

            if (alphaVal > metalFractionTol)
            {
                inFilm = true;
            }
            else if (alphaVal < metalFractionTol)
            {
                inFilm = false;
            }
        }
        bool depositInMetal = inFilm;
        const scalar gasCoeff = gasAbsorptionCoeff_.value();

        if (metalFractions)
        {
            const scalar alphaVal = Foam::min
            (
                Foam::max((*metalFractions)[cellI], scalar(0)),
                scalar(1)
            );

            if (alphaVal > metalFractionTol)
            {
                depositInMetal = true;
            }
            else
            {
                depositInMetal = false;
            }
        }
        if (!isInBeam(c, axialHalfLength))
        {
            continue;
        }

        ++metrics.cellsInBeam;
        metrics.totalBeamVolume += cellVolumes[cellI];
        if (depositInMetal)
        {
            ++metrics.cellsInFilm;
        }
        else
        {
            ++metrics.cellsInGas;

            if (gasCoeff <= VSMALL)
            {
                continue;
            }
        }
        vector r = c - focus_;
        scalar z = (r & direction_);
        r -= z*direction_;
        const scalar R = mag(r);

        const scalar spatialTerm = calculateGaussianIntensity(R);

        const scalar cellAbsorptionCoeff =
            depositInMetal ? absorptionCoeff_.value() : gasCoeff;

        const labelList& pointLabels = mesh_.cellPoints(cellI);

        if (pointLabels.empty())
        {
            continue;
        }

        scalar sMin = GREAT;
        scalar sMax = -GREAT;

        forAll(pointLabels, pointi)
        {
            const label ptI = pointLabels[pointi];
            const scalar s = (meshPoints[ptI] - entryPoint) & directionUnit;
            sMin = Foam::min(sMin, s);
            sMax = Foam::max(sMax, s);
        }

        scalar sIn = sMin;
        scalar sOut = sMax;

        if (sOut < sIn)
        {
            const scalar tmp = sIn;
            sIn = sOut;
            sOut = tmp;
        }
        sIn = Foam::max(sIn, scalar(0));
        sOut = Foam::max(sOut, sIn);
        const bool cellUsesFilmInterval = inFilm && haveFilmInterval;

        if (cellUsesFilmInterval)
        {
            sIn = Foam::max(sIn, filmIntervalStart);
            sOut = Foam::min(sOut, filmIntervalEnd);

            if (sOut <= sIn)
            {
                continue;
            }
        }        
        // reflectivity_ is defined at the film interface, so the transmission
        // derived from it attenuates the beam prior to entering either phase.

        scalar transmissionFactor = interfaceTransmission;

        if (!applyInterfaceTransmissionEverywhere_ && !inFilm)
        {
            transmissionFactor = 1.0;
        }

        const scalar baseIntensity =
              peakIntensity_.value()
            * temporalAverage
            * spatialTerm
            * transmissionFactor;

        // ========== PLASMA SHIELDING ==========
        // Check for plasma shielding (reduces laser absorption)
        scalar shieldingFactor = 1.0;
        if (mesh_.foundObject<volScalarField>("plasmaShielding"))
        {
            const volScalarField& plasmaShield =
                mesh_.lookupObject<volScalarField>("plasmaShielding");

            // shieldingFactor ranges from 0 (full shielding) to 1 (no shielding)
            shieldingFactor = Foam::max(1.0 - plasmaShield[cellI], scalar(0.0));
        }

        // Apply plasma shielding to laser intensity
        const scalar effectiveIntensity = baseIntensity * shieldingFactor;
        // ========== END PLASMA SHIELDING ==========

        const scalar deltaS = Foam::max(sOut - sIn, VSMALL);
        scalar Ein = effectiveIntensity;
        scalar Eout = effectiveIntensity;

        scalar sInForExponent = sIn;
        scalar sOutForExponent = sOut;

        if (inFilm)
        {
            const scalar depthIn = Foam::max(scalar(0), sIn - filmEntryOffset);
            const scalar depthOut = Foam::max(depthIn, sOut - filmEntryOffset);

            sInForExponent = depthIn;
            sOutForExponent = depthOut;
        }

        if (cellAbsorptionCoeff > VSMALL)
        {
            Ein *= exp(-cellAbsorptionCoeff*sInForExponent);
            Eout *= exp(-cellAbsorptionCoeff*sOutForExponent);
        }
        scalar qVol = (Ein - Eout)/deltaS;
        qVol = Foam::max(qVol, scalar(0));

        if (std::isfinite(qVol) && qVol > 0)
        {
            scalar limitedValue = qVol;

            if (metrics.limitSource && limitedValue > metrics.maxSourceCap)
            {
                limitedValue = metrics.maxSourceCap;
                ++metrics.limitedCells;
            }

            limitedValue = Foam::max(limitedValue, scalar(0));

            source[cellI] = limitedValue;

            const scalar cellPower = source[cellI] * cellVolumes[cellI];

            metrics.maxSourceValue =
                max(metrics.maxSourceValue, source[cellI]);
            metrics.totalSourceIntegral += cellPower;

            if (depositInMetal)
            {
                metrics.totalFilmSourceIntegral += cellPower;
            }
            else
            {
                metrics.totalGasSourceIntegral  += cellPower;
            }
        }
    }
    // After beam weighting loop
    if (verbose && Pstream::master())
    {
        Info<< "BEAM PROFILE VALIDATION:" << nl
            << "  Total mesh cells: " << cellCentres.size() << nl
            << "  Cells in search box: " << cellsInSearchBox << nl
            << "  Cells in beam: " << metrics.cellsInBeam << nl
            << "  Peak source: " << metrics.maxSourceValue << " W/m³" << nl
            << "  Integrated power: " << metrics.totalSourceIntegral << " W" << nl
            << "  Average source: " << metrics.totalSourceIntegral/metrics.totalBeamVolume
            << " W/m³" << nl
            << "  Peak/Average ratio: " << metrics.maxSourceValue * metrics.totalBeamVolume
               / metrics.totalSourceIntegral << nl;
    }
    reduce(metrics.cellsInBeam, sumOp<label>());
    reduce(metrics.cellsInFilm, sumOp<label>());
    reduce(metrics.cellsInGas,  sumOp<label>());
    reduce(metrics.maxSourceValue, maxOp<scalar>());
    reduce(metrics.totalSourceIntegral, sumOp<scalar>());
    reduce(metrics.totalBeamVolume, sumOp<scalar>());
    reduce(metrics.totalFilmSourceIntegral, sumOp<scalar>());
    reduce(metrics.totalGasSourceIntegral, sumOp<scalar>());
    reduce(metrics.limitedCells, sumOp<label>());

    return metrics;
    
}

//------------------------------------------------------------------------------
void femtosecondLaserModel::emitDiagnostics
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
) const
{
    if (!verbose)
    {
        return;
    }

    if (metrics.limitSource && metrics.limitedCells > 0)
    {
        Info<< "Laser source limited to " << metrics.maxSourceCap
            << " W/m^3 in " << metrics.limitedCells << " cells" << endl;
    }

    const scalar avgIntensityInBeam =
        (metrics.totalBeamVolume > 0)
      ? metrics.totalSourceIntegral/metrics.totalBeamVolume
      : 0.0;

    Info<< "LASER DIAGNOSTICS:" << nl
        << "  Input peak intensity: " << peakIntensity_.value() << " W/m^2" << nl
        << "  Average volumetric intensity in beam: " << avgIntensityInBeam << " W/m^3" << nl
        << "  Absorption coefficient (film): " << absorptionCoeff_.value() << " 1/m" << nl
        << "  Absorption coefficient (gas):  " << gasAbsorptionCoeff_.value() << " 1/m" << nl
        << "  Spot radius: " << spotSize_.value()/2.0*1e6 << " µm" << nl
        << "  Beam area: "
        << constant::mathematical::pi*sqr(spotSize_.value()/2.0)*1e12
        << " µm^2" << endl;

    if (depositionActive && (timeIndex % 10 == 0))
    {
        const scalar dt = dtDim.value();

        Info<< "LASER ENERGY DEPOSITION:" << nl
            << "  Time: " << currentTime*1e12 << " ps" << nl
            << "  Temporal factor: " << envelope.temporalAverage << nl
            << "  Cells in beam: " << metrics.cellsInBeam << nl
            << "  Cells in metal film: " << metrics.cellsInFilm << nl
            << "  Cells in gas: " << metrics.cellsInGas << nl
            << "  Max intensity: " << metrics.maxSourceValue << " W/m^3" << nl
            << "  Total power: " << metrics.totalSourceIntegral << " W" << nl
            << "    Film power: " << metrics.totalFilmSourceIntegral << " W" << nl
            << "    Gas power:  " << metrics.totalGasSourceIntegral  << " W" << nl
            << "  dt: " << dt << " s" << nl
            << "  Incident energy this step: " << incidentEnergyThisStep << " J" << nl
            << "  Absorbed energy this step: " << energyThisStep.value() << " J" << nl
            << "    Film energy this step: " << filmEnergyThisStep << " J" << nl
            << "    Gas energy this step:  " << gasEnergyThisStep  << " J" << nl
            << "  Cumulative incident energy: " << cumulativeIncidentEnergy_ << " J" << nl
            << "  Cumulative absorbed energy: " << cumulativeEnergy_ << " J" << nl
            << "    Cumulative film: " << cumulativeFilmEnergy_ << " J" << nl
            << "    Cumulative gas:  " << cumulativeGasEnergy_  << " J" << endl;

        if (metrics.cellsInBeam == 0)
        {
            WarningInFunction
                << "No cells in laser beam!"
                << " Focus: " << focus_
                << " Bounds: " << mesh_.bounds() << endl;
        }
        if (metrics.cellsInFilm == 0 && metrics.cellsInBeam > 0)
        {
            WarningInFunction
                << "Beam hits " << metrics.cellsInBeam
                << " cells but none in metal film region" << endl;
        }
        else if (metrics.cellsInGas > 0 && gasAbsorptionCoeff_.value() <= VSMALL)
        {
            Info<< "    Gas absorption disabled (gasAbsorptionCoeff ≈ 0);"
                << " no direct heating deposited off-film." << endl;
        }
    }
}

//------------------------------------------------------------------------------
void femtosecondLaserModel::updateEnergyTracking
(
    const volScalarField& source,
    const dimensionedScalar& dtDim,
    const EnvelopeResult& envelope,
    const SpatialMetrics& metrics,
    const scalar overlapStart,
    const scalar currentTime,
    const bool depositionActive,
    const label timeIndex
) const
{
    const scalar dt = dtDim.value();

    const dimensionedScalar energyThisStep =
        fvc::domainIntegrate(source * dtDim);
    const scalar stepEnergy = energyThisStep.value();
    const scalar incidentEnergyThisStep = envelope.expectedEnergy;
    cumulativeEnergy_ += stepEnergy;
    cumulativeIncidentEnergy_ += incidentEnergyThisStep;

    if (!continuousLaser_)
    {
        const scalar depositableFraction = depositableEnergyFraction();
        const scalar coverageFraction = beamCoverageFraction();
        const scalar depositableExpectedThisStep =
            incidentEnergyThisStep*depositableFraction*coverageFraction;
        const scalar rawStepEnergy = metrics.rawSourceIntegral*dt;

        if (depositionActive || trackingPulse_)
        {
            if (!trackingPulse_)
            {
                trackingPulse_ = true;
                pulseCompleted_ = false;
                pulseEnergyAccumulator_ = 0.0;
                pulseExpectedAccumulator_ = 0.0;
                pulseDepositableExpectedAccumulator_ = 0.0;
                pulseRawEnergyAccumulator_ = 0.0;
                currentPulseStartTime_ = overlapStart;
            }

            pulseEnergyAccumulator_ += stepEnergy;
            pulseExpectedAccumulator_ += incidentEnergyThisStep;
            pulseDepositableExpectedAccumulator_ += depositableExpectedThisStep;
            pulseRawEnergyAccumulator_ += rawStepEnergy;
        }
    }

    const scalar filmEnergyThisStep = metrics.totalFilmSourceIntegral * dt;
    const scalar gasEnergyThisStep  = metrics.totalGasSourceIntegral  * dt;

    cumulativeFilmEnergy_ += filmEnergyThisStep;
    cumulativeGasEnergy_  += gasEnergyThisStep;

    if (!checkEnergyConservation())
    {
        WarningInFunction
            << "Energy check failed at step " << timeIndex
            << " (E=" << energyThisStep.value() << " J)" << endl;
    }

    emitDiagnostics
    (
        envelope,
        metrics,
        dtDim,
        energyThisStep,
        incidentEnergyThisStep,
        filmEnergyThisStep,
        gasEnergyThisStep,
        currentTime,
        timeIndex,
        depositionActive
    );
}
//------------------------------------------------------------------------------
// source compute (temporal envelope + spatial weighting)
//------------------------------------------------------------------------------
void femtosecondLaserModel::calculateSource() const
{
    if (sourceValid_)
    {
        return;
    }

    activeThisStep_ = false;
    volScalarField& source = resetSourceField();

    const scalar t = mesh_.time().value();
    const dimensionedScalar dtDim = mesh_.time().deltaT();
    const scalar dt = dtDim.value();
    const scalar tStart = t - dt;
    const label timeIndex = mesh_.time().timeIndex();
    const bool master = Pstream::master();
    // ENHANCED DEBUG OUTPUT
    if (verbose && master && timeIndex % 10 == 0)
    {
         Info<< "=== LASER SOURCE CALCULATION ENTERED ===" << nl
            << "Time: " << mesh_.time().value() << " s" << nl
            << "Focus: " << focus_ << nl
            << "Mesh bounds: " << mesh_.bounds() << nl;
                   
        Info<< "===== LASER DEBUG =====" << nl
            << "Time: " << t << " s (" << t*1e12 << " ps)" << nl
            << "Time step: " << dt << " s (" << dt*1e12 << " ps)" << nl
            << "Laser window: [" << laserStartTime_ << ", " << laserEndTime_ << "] s" << nl
            << "Focus: " << focus_ << nl
            << "  Focus Y: " << focus_.y()*1e6 << " µm" << nl
            << "Film bounds: [" << filmYMin_*1e6 << ", " << filmYMax_*1e6 << "] µm" << nl
            << "Peak intensity: " << peakIntensity_.value() << " W/m^2" << nl
            << "Pulse energy: " << pulseEnergy_.value() << " J" << nl
            << "Spot size: " << spotSize_.value()*1e6 << " µm diameter" << nl
            << "Direction: " << direction_ << nl
            << "Absorption coeff: " << absorptionCoeff_.value() << " 1/m" << endl;
    }

    if (timeIndex < lastTimeIndex_)
    {
        cumulativeEnergy_ = 0.0;
        cumulativeIncidentEnergy_ = 0.0;
        cumulativeFilmEnergy_ = 0.0;
        cumulativeGasEnergy_  = 0.0;
        trackingPulse_ = false;
        pulseCompleted_ = false;
        pulseEnergyAccumulator_ = 0.0;
        pulseExpectedAccumulator_ = 0.0;
        pulseDepositableExpectedAccumulator_ = 0.0;
        pulseRawEnergyAccumulator_ = 0.0;
        currentPulseStartTime_ = 0.0;
        pulseCounter_ = 0;
    }
    lastTimeIndex_ = timeIndex;

    const bool singlePulse = (!continuousLaser_ && pulseFrequency_ <= SMALL);
    const bool multiPulse  = (!continuousLaser_ && pulseFrequency_ > SMALL);

    if (singlePulse && pulseCompleted_)
    {
        sourceValid_ = true;
        return;
    }

    const scalar overlapStart = max(tStart, laserStartTime_);
    const scalar overlapEnd   = min(t, laserEndTime_);
    const scalar overlapDuration = Foam::max(0.0, overlapEnd - overlapStart);
    const scalar comparisonLaserEnd = laserEndTime_;

    if (verbose && master && timeIndex % 10 == 0)
    {
        Info<< "Time overlap check:" << nl
            << "  Current time window: [" << tStart*1e12 << ", " << t*1e12 << "] ps" << nl
            << "  Laser active window: [" << laserStartTime_*1e12 << ", " << laserEndTime_*1e12 << "] ps" << nl
            << "  Overlap window: [" << overlapStart*1e12 << ", " << overlapEnd*1e12 << "] ps" << nl
            << "  Overlap duration: " << overlapDuration*1e12 << " ps" << endl;
    }

    if (overlapDuration <= VSMALL)
    {
        const scalar tol = 1e-12;
        const bool beforePulse = t <= laserStartTime_ + tol;
        const bool afterPulse = tStart >= comparisonLaserEnd - tol;

        if (beforePulse)
        {
            if (verbose && master && timeIndex % 10 == 0)
            {
                Info<< "No time overlap - before laser start" << endl;
            }
            sourceValid_ = true;
            return;
        }

        if (afterPulse)
        {
            if (verbose && master && timeIndex % 10 == 0)
            {
                Info<< "No time overlap - laser pulse finished" << endl;
            }
            if (!pulseCompleted_)
            {
                finalizePulseEnergyCheck("post-pulse inactive window", t);
                pulseCompleted_ = true;
            }
            sourceValid_ = true;
            return;
        }

        if (verbose && master && timeIndex % 10 == 0)
        {
            Info<< "No time overlap - inside laser window without overlap" << endl;
        }
        sourceValid_ = true;
        return;
    }

    EnvelopeResult envelope =
        evaluateTemporalEnvelope(overlapStart, overlapStart + overlapDuration, dt);

    if (verbose && master && timeIndex % 10 == 0)
    {
        Info<< "Temporal envelope:" << nl
            << "  Active: " << (envelope.active ? "YES" : "NO") << nl
            << "  Temporal integral: " << envelope.temporalIntegral << " s" << nl
            << "  Temporal average: " << envelope.temporalAverage << nl
            << "  Expected incident energy: " << envelope.expectedEnergy
            << " J" << endl; // coverage ignores tree-box clipping near bounds
    }

    scalar sigma = 0.0;
    scalar pulseCenter = laserStartTime_;

    if (singlePulse)
    {
        sigma = pulseWidth_.value()/(2.0*sqrt(2.0*log(2.0)));

        if (pulseCenterTime_ > -0.5*GREAT)
        {
            pulseCenter = pulseCenterTime_;
        }
        else
        {
            const scalar windowWidth = max(laserEndTime_ - laserStartTime_, SMALL);

            if (sigma <= VSMALL)
            {
                pulseCenter = laserStartTime_ + 0.5*windowWidth;
            }
            else
            {
                const scalar maxCenterOffset = 3.0*sigma;
                const scalar halfWindow = 0.5*windowWidth;
                const scalar centerOffset = min(halfWindow, maxCenterOffset);
                pulseCenter = laserStartTime_ + centerOffset;
            }
        }
    }

    const scalar pulseToleranceRadius = (singlePulse && sigma > VSMALL)
        ? 5.0*sigma
        : 0.0;
    const scalar toleranceStart = pulseCenter - pulseToleranceRadius;
    const scalar toleranceEnd = pulseCenter + pulseToleranceRadius;

    if (!envelope.active || envelope.temporalAverage <= VSMALL)
    {
        bool outsidePulseWindow = false;

        if (singlePulse)
        {
            if (sigma <= VSMALL)
            {
                outsidePulseWindow = true;
            }
            else
            {
                const bool rawOutside =
                    (t <= toleranceStart || tStart >= toleranceEnd);

                outsidePulseWindow =
                    rawOutside || overlapStart >= toleranceEnd;
            }
        }

        const bool haveOverlap = overlapDuration > VSMALL;

        bool withinActiveWindow =
            haveOverlap
         && t >= laserStartTime_
         && t <= comparisonLaserEnd;
        if (withinActiveWindow && multiPulse)
        {
            bool overlapsPulseWindow = false;

            const scalar period = 1.0/pulseFrequency_;
            const scalar onTime = max(SMALL, pulseDutyCycle_*period);

            if (onTime > VSMALL)
            {
                const scalar localStart = Foam::max(0.0, overlapStart - laserStartTime_);
                const scalar localEnd   = Foam::max(localStart, overlapEnd - laserStartTime_);

                label periodIndex = 0;

                if (localStart > SMALL)
                {
                    periodIndex = static_cast<label>(std::floor(localStart/period));
                }

                for
                (
                    scalar windowStart = periodIndex*period;
                    windowStart < localEnd + period;
                    windowStart += period
                )
                {
                    const scalar windowEnd = windowStart + onTime;

                    if (windowEnd <= localStart)
                    {
                        continue;
                    }

                    if (windowStart >= localEnd)
                    {
                        break;
                    }

                    overlapsPulseWindow = true;
                    break;
                }
            }

            if (!overlapsPulseWindow)
            {
                withinActiveWindow = false;
            }
        }
        if (singlePulse && outsidePulseWindow)
        {
            withinActiveWindow = false;
            if (!pulseCompleted_)
            {
                WarningInFunction
                    << "Time window [" << tStart*1e12 << ", " << t*1e12
                    << "] ps lies outside the ±5σ Gaussian pulse envelope"
                    << " centred at " << pulseCenter*1e12 << " ps. "
                    << "Laser deposition stops even though laserEndTime = "
                    << laserEndTime_*1e12 << " ps. "
                    << "Increase pulseCenterTime or pulseFrequency if energy"
                    << " is required later in the simulation." << endl;
            }
            pulseCompleted_ = true;
        }


        if
        (
            withinActiveWindow
         && !(singlePulse && t > comparisonLaserEnd + 1e-12)
        )
        {
           WarningInFunction
                << "Temporal envelope is inactive while current time "
                << "slice lies inside the configured laser window. "
                << "No laser energy will be deposited." << endl;
        }

        if (verbose && master && timeIndex % 10 == 0)
        {
            Info<< "Temporal envelope inactive - no laser heating" << endl;
        }
        finalizePulseEnergyCheck("pulse window complete", t);
        sourceValid_ = true;
        return;
    }

    if (verbose && master && timeIndex % 10 == 0)
    {
        Info<< "Applying spatial weighting..." << nl
            << "  Temporal average factor: " << envelope.temporalAverage << endl;

        const boundBox& bounds = mesh_.bounds();
        const bool focusInMesh = bounds.contains(focus_);
        Info<< "Geometry checks:" << nl
            << "  Mesh bounds: " << bounds << nl
            << "  Focus position: " << focus_ << nl
            << "  Focus in mesh: " << (focusInMesh ? "YES" : "NO") << nl;

        if (!focusInMesh)
        {
            WarningInFunction
                << "Laser focus is outside mesh bounds!" << nl
                << "This will result in zero heating." << endl;
        }

        const bool focusInFilm = (focus_.y() >= filmYMin_ && focus_.y() <= filmYMax_);
        Info<< "  Focus Y: " << focus_.y()*1e6 << " µm" << nl
            << "  Film Y bounds: [" << filmYMin_*1e6 << ", " << filmYMax_*1e6 << "] µm" << nl
            << "  Focus in film: " << (focusInFilm ? "YES" : "NO") << endl;

        if (!focusInFilm)
        {
            WarningInFunction
                << "Laser focus Y-position is outside film region!" << nl
                << "This may result in reduced or zero film heating." << endl;
        }
    }

    SpatialMetrics metrics = applySpatialWeighting(source, envelope.temporalAverage);

    const scalar dtValue = dtDim.value();

    if (envelope.expectedEnergy <= VSMALL || dtValue <= VSMALL)
    {
        if (metrics.totalSourceIntegral > VSMALL)
        {
            source = dimensionedScalar("zero", source.dimensions(), 0.0);
            metrics.totalSourceIntegral = 0.0;
            metrics.totalFilmSourceIntegral = 0.0;
            metrics.totalGasSourceIntegral = 0.0;
            metrics.maxSourceValue = 0.0;
        }
    }
    else
    {
        const scalar incidentPower = envelope.expectedEnergy/max(dtValue, VSMALL);
        const scalar depositableFraction = depositableEnergyFraction();
        const scalar coverageFraction = beamCoverageFraction();
        const scalar desiredPower =
            incidentPower*depositableFraction*coverageFraction;
        const scalar actualPower = metrics.totalSourceIntegral;
        metrics.rawSourceIntegral = actualPower;

        if (actualPower > VSMALL)
        {
            scalar scale = desiredPower/actualPower;

            if (metrics.limitSource && metrics.maxSourceValue > VSMALL)
            {
                const scalar maxScale = metrics.maxSourceCap/metrics.maxSourceValue;

                if (scale > maxScale)
                {
                    scale = maxScale;
                }
            }

            if (mag(scale - 1.0) > 1e-6)
            {
                source *= scale;
                metrics.totalSourceIntegral      *= scale;
                metrics.totalFilmSourceIntegral *= scale;
                metrics.totalGasSourceIntegral  *= scale;
                metrics.maxSourceValue          *= scale;
            }
        }
        else
        {
            WarningInFunction
                << "Temporal envelope expects " << envelope.expectedEnergy
                << " J this step, but spatial integration returned zero"
                << " power.  No laser energy will be deposited." << endl;
        }
    }

    if (verbose && master && timeIndex % 10 == 0)
    {
        Info<< "Spatial processing results:" << nl
            << "  Cells in beam: " << metrics.cellsInBeam << nl
            << "  Cells in film: " << metrics.cellsInFilm << nl
            << "  Cells in gas: " << metrics.cellsInGas << nl
            << "  Max source value: " << metrics.maxSourceValue << " W/m^3" << nl
            << "  Total power: " << metrics.totalSourceIntegral << " W" << nl
            << "  Film power: " << metrics.totalFilmSourceIntegral << " W" << nl
            << "  Limited cells: " << metrics.limitedCells << endl;
    }

    if (metrics.cellsInBeam == 0)
    {
        WarningInFunction
            << "No cells found in laser beam path at time " << t*1e12 << " ps!" << nl
            << "Check focus position, spot size, and mesh resolution." << nl
            << "DIAGNOSTIC INFO:" << nl
            << "  Focus position: " << focus_ << nl
            << "  Beam direction: " << direction_ << nl
            << "  Spot size: " << spotSize_.value() << " m" << nl
            << "  Film Y bounds: [" << filmYMin_ << ", " << filmYMax_ << "]" << nl
            << "  Mesh bounds: [" << mesh_.bounds().min() << ", " << mesh_.bounds().max() << "]" << nl
            << "  Total mesh cells: " << mesh_.nCells() << nl
            << endl;
    }
    else if (metrics.cellsInFilm == 0 && metrics.cellsInBeam > 0)
    {
        WarningInFunction
            << "Beam intersects " << metrics.cellsInBeam
            << " cells but none are in film region!" << nl
            << "Check film bounds and focus Y-position." << endl;
    }
    else if (metrics.maxSourceValue <= VSMALL && metrics.cellsInBeam > 0)
    {
        WarningInFunction
            << "Beam intersects " << metrics.cellsInBeam
            << " cells but produces zero heating!" << nl
            << "Check absorption coefficients and temporal averaging." << endl;
    }

    activeThisStep_ = (metrics.maxSourceValue > VSMALL);

    if (verbose && master && timeIndex % 10 == 0)
    {
        Info<< "Laser active this step: " << (activeThisStep_ ? "YES" : "NO") << nl;
        if (activeThisStep_)
        {
            Info<< "  SUCCESS: Laser is depositing energy!" << nl;
        }
        Info<< "=====================" << endl;
    }
    updateEnergyTracking
    (
        source,
        dtDim,
        envelope,
        metrics,
        overlapStart,
        t,
        activeThisStep_,
        timeIndex
    );

    sourceValid_ = true;
}

//------------------------------------------------------------------------------
tmp<volScalarField> femtosecondLaserModel::source() const
{
    calculateSource();
    return tSource_;
}
//------------------------------------------------------------------------------
bool femtosecondLaserModel::activeThisStep() const
{
    calculateSource();
    return activeThisStep_;
}

//------------------------------------------------------------------------------
bool femtosecondLaserModel::valid() const
{
    validateParameters();
    return checkPhysicalBounds();
}

//------------------------------------------------------------------------------
void femtosecondLaserModel::write() const
{
    const bool master = Pstream::master();
    if (verbose && master)
    {
        Info<< "Femtosecond laser model status:" << nl
            << "  Mode: " << (continuousLaser_ ? "Continuous" : "Pulsed") << nl
            << "  Peak intensity: " << peakIntensity_.value() << " W/m^2" << nl
            << "  Pulse width: " << pulseWidth_.value() << " s" << nl
            << "  Wavelength: " << wavelength_.value() << " m" << nl
            << "  Spot size: " << spotSize_.value() << " m" << nl
            << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
            << "  Max volumetric source cap: " << maxVolumetricSource_.value()
            << " W/m^3" << nl;

        if (tSource_.valid())
        {
            const dimensionedScalar currentMax = max(tSource_());
            Info<< "  Current max volumetric source: "
                << currentMax.value() << " W/m^3" << nl;
        }

        Info<< "  Absorption coefficient: " << absorptionCoeff_.value() << " 1/m" << nl
            << "  Gas absorption coefficient: " << gasAbsorptionCoeff_.value() << " 1/m" << nl
            << "  Reflectivity: " << reflectivity_ << nl
            << "  Transmission: " << (transmission_ >= 0
                                      ? transmission_ : (1.0 - reflectivity_)) << nl
            << "  Incidence angle: " << incidenceAngle_ << " rad" << nl
            << "  Focus: " << focus_ << nl
            << "  Direction: " << direction_ << nl
            << "  Active time: " << laserStartTime_ << " to "
            << laserEndTime_ << " s" << endl;

        if (tSource_.valid())
        {
            const dimensionedScalar maxI = max(tSource_());
            const dimensionedScalar e    =
                fvc::domainIntegrate(tSource_()*mesh_.time().deltaT());

            Info<< "Source stats: max=" << maxI.value() << " W/m^3, "
                << "E(step)=" << e.value() << " J" << nl
                << "  Cumulative incident: " << cumulativeIncidentEnergy_ << " J" << nl
                << "  Cumulative absorbed: " << cumulativeEnergy_ << " J" << nl
                << "    Film cumulative: " << cumulativeFilmEnergy_ << " J" << nl
                << "    Gas cumulative:  " << cumulativeGasEnergy_  << " J" << endl;
        }
    }
}

} // namespace Foam
