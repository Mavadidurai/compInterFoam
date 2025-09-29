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
    filmYMin_(0.0),
    filmYMax_(0.0),
    tSource_(),
    sourceValid_(false),
    cumulativeEnergy_(0.0),
    cumulativeFilmEnergy_(0.0),
    cumulativeGasEnergy_(0.0),
    lastTimeIndex_(mesh.time().timeIndex()),
    pulseEnergyAccumulator_(0.0),
    pulseExpectedAccumulator_(0.0),
    currentPulseStartTime_(0.0),
    pulseCounter_(0),
    trackingPulse_(false),
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
            << "  Max volumetric source: " << maxVolumetricSource_.value()
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
bool femtosecondLaserModel::validateParameters() const
{
    bool ok = true;

    ok = ok && (peakIntensity_.dimensions() == dimPower/dimArea);
    ok = ok && (pulseWidth_.dimensions()    == dimTime);
    ok = ok && (wavelength_.dimensions()    == dimLength);
    ok = ok && (absorptionCoeff_.dimensions()== dimless/dimLength);
    ok = ok && (gasAbsorptionCoeff_.dimensions()== dimless/dimLength);    
    ok = ok && (spotSize_.dimensions()      == dimLength);
    ok = ok && (pulseEnergy_.dimensions()   == dimEnergy);
    ok = ok && (maxVolumetricSource_.dimensions() == dimPower/dimVolume);

    ok = ok && (peakIntensity_.value() > 0);
    ok = ok && (pulseWidth_.value()    > 0);
    ok = ok && (wavelength_.value()    > 0);
    ok = ok && (spotSize_.value()      > 0);
    ok = ok && (pulseEnergy_.value()   > 0);
    ok = ok && (gasAbsorptionCoeff_.value() >= 0);
    ok = ok && (maxVolumetricSource_.value() >= 0);
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
    return ok;
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
        currentPulseStartTime_ = 0.0;
        return;
    }

    const scalar expected = pulseExpectedAccumulator_;
    const scalar configuredDepositable =
        pulseEnergy_.value()*depositableEnergyFraction();
    const scalar reference =
        max(max(expected, configuredDepositable), scalar(1e-12));
    const scalar tolerance =
        max(pulseEnergyToleranceAbs_, pulseEnergyToleranceRel_*reference);
    const scalar diffExpected = mag(pulseEnergyAccumulator_ - expected);
    const scalar diffConfigured =
        mag(pulseEnergyAccumulator_ - configuredDepositable);
    const scalar expectedMismatch =
        mag(expected - configuredDepositable);
    const scalar configuredResidual =
        max(diffConfigured - expectedMismatch, scalar(0));    
    const scalar maxDeviation = max(diffExpected, diffConfigured);

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
            << "  Expected (integrated):      " << expected << " J" << nl
            << "  Configured depositable:     " << configuredDepositable << " J"
            << nl
            << "  Diff vs expected:          " << diffExpected << " J" << nl
            << "  Diff vs configured depositable:   " << diffConfigured
            << " J" << nl
            << "  Expected-configured depositable diff: "
            << expectedMismatch << " J" << nl
            << "  Configured depositable residual:  " << configuredResidual
            << " J" << nl
            << "  Max deviation:        " << maxDeviation << " J" << nl
            << "  Tolerance(abs):       " << pulseEnergyToleranceAbs_ << " J" << nl
            << "  Tolerance(rel*ref):   "
            << pulseEnergyToleranceRel_*reference << " J" << nl
            << "  Applied tolerance:    " << tolerance << " J" << endl;
    }

    const bool warnConfigured = configuredResidual > tolerance;
    const bool warnExpected = diffExpected > tolerance;

    if (warnConfigured)
    {
        WarningInFunction
            << "Laser pulse energy mismatch after pulse " << pulseCounter_
            << ": deposited " << pulseEnergyAccumulator_
            << " J vs configured depositable energy " << configuredDepositable
            << " J (reflective/transmissive losses excluded). Residual mismatch"
            << " " << configuredResidual << " J exceeds tolerance "
            << tolerance << " J" << endl;
    }

    if (warnExpected)
    {
        WarningInFunction
            << "Laser pulse energy mismatch after pulse " << pulseCounter_
            << ": deposited " << pulseEnergyAccumulator_ << " J vs integrated"
            << " expectation " << expected << " J (difference "
            << diffExpected << " J exceeds tolerance " << tolerance << " J)"
            << endl;
    }

    if (!warnConfigured && !warnExpected && verbose && master)
    {
        if (diffConfigured > VSMALL)
        {
            Info<< "  Note: deposited energy differs from configured depositable"
                << " energy by " << diffConfigured
                << " J (within tolerance)." << endl;
        }

        if (diffExpected > VSMALL)
        {
            Info<< "  Note: deposited energy differs from integrated expectation"
                << " by " << diffExpected << " J (within tolerance)." << endl;
        }

        if (expectedMismatch > VSMALL && expectedMismatch <= tolerance)
        {
            Info<< "  Note: integrated expectation differs from configured"
                << " depositable energy by " << expectedMismatch
                << " J (within tolerance)." << endl;
        }
    }

    trackingPulse_ = false;
    pulseEnergyAccumulator_ = 0.0;
    pulseExpectedAccumulator_ = 0.0;
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
    if (sigma <= VSMALL)
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
    else if (incidenceAngle_ > VSMALL)
    {
        const scalar n1 = 1.0;
        const scalar sqrtR = sqrt(max(reflectivity, scalar(0)));
        const scalar n2 = (1.0 + sqrtR)/max(VSMALL, (1.0 - sqrtR));
        const scalar sinThetaT = n1/n2 * sin(incidenceAngle_);

        if (mag(sinThetaT) < 1.0)
        {
            const scalar cosTheta  = cos(incidenceAngle_);
            const scalar cosThetaT = sqrt(1.0 - sqr(sinThetaT));
            const scalar Rs = sqr((n1*cosTheta - n2*cosThetaT)
                                /(n1*cosTheta + n2*cosThetaT));
            const scalar Rp = sqr((n1*cosThetaT - n2*cosTheta)
                                /(n1*cosThetaT + n2*cosTheta));
            transmissionFactor = 1.0 - 0.5*(Rs + Rp);
        }
        else
        {
            transmissionFactor = 0.0;
        }
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

    return coverage;
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
    const scalar depositableFraction = depositableEnergyFraction();
    const scalar coverageFraction = beamCoverageFraction();    
    if (continuousLaser_)
    {
        result.temporalIntegral = overlapEnd - overlapStart;
    }
    else if (pulseFrequency_ > SMALL)
    {

        const scalar period   = 1.0/pulseFrequency_;
        const scalar onTime   = max(SMALL, pulseDutyCycle_*period);
        const scalar sigma    = pulseWidth_.value()/(2.0*sqrt(2.0*log(2.0)));
        const scalar localStart = overlapStart - laserStartTime_;
        const scalar localEnd   = overlapEnd   - laserStartTime_;

        label periodIndex = 0;
        if (localStart > SMALL)
        {
            periodIndex = static_cast<label>(std::floor(localStart/period));
        }
        // Normalised single-pulse integral used for energy expectations
        const scalar fullPulseIntegral =
            max
            (
                VSMALL,
                2.0*sigma*sqrt(constant::mathematical::pi/2.0)
              * std::erf(0.5*onTime/(sqrt(2.0)*sigma))
            );

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
            const scalar integral =
                gaussianWindowIntegral(clipStart, clipEnd, center, sigma);

            if (integral > VSMALL)
            {
                result.temporalIntegral += integral;
            }
        }

        if (result.temporalIntegral > VSMALL)
        {
            result.expectedEnergy =
                pulseEnergy_.value()
              * depositableFraction
              * coverageFraction    
              * result.temporalIntegral/max(fullPulseIntegral, VSMALL);
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
        // configured start time while preserving backward compatibility when the
        // time window is narrow.

        const scalar sigma =
            pulseWidth_.value()/(2.0*sqrt(2.0*log(2.0)));  // std dev from FWHM
        const scalar windowWidth = max(laserEndTime_ - laserStartTime_, SMALL);

        // Limit how far the peak can move away from the trigger time.  Using
        // ±3σ retains >99% of the Gaussian energy while avoiding large delays
        // when laserEndTime is set long after the pulse.
        const scalar maxCenterOffset = 3.0*sigma;
        const scalar halfWindow = 0.5*windowWidth;
        const scalar centerOffset = min(halfWindow, maxCenterOffset);
        const scalar center = laserStartTime_ + centerOffset;

        result.temporalIntegral =
            gaussianWindowIntegral(overlapStart, overlapEnd, center, sigma);

        if (result.temporalIntegral > VSMALL)
        {
            const scalar fullIntegral = sigma*sqrt(2.0*constant::mathematical::pi);
            result.expectedEnergy =
                pulseEnergy_.value()
              * depositableFraction
              * coverageFraction               
              * result.temporalIntegral/max(fullIntegral, VSMALL);
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

    const scalar beamRadius = spotSize_.value()/2.0;
    const scalar radialHalfWidth = 3.0*beamRadius; // ~3-sigma laterally

    const scalar filmHalfThickness =
        0.5*Foam::max(filmYMax_ - filmYMin_, scalar(0));

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

    const pointField& cellCentres = mesh_.C();
    const scalarField& cellVolumes = mesh_.V();
    const pointField& meshPoints = mesh_.points();

    vector directionUnit(direction_);
    const scalar directionMag = mag(directionUnit);

    if (directionMag > VSMALL)
    {
        directionUnit /= directionMag;
    }
    else
    {
        directionUnit = vector::zero;
    }

    point entryPoint = focus_;

    if (directionMag > VSMALL)
    {
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
    }
    scalar filmIntervalStart = 0.0;
    scalar filmIntervalEnd = 0.0;
    bool haveFilmInterval = false;

    if (directionMag > VSMALL)
    {
        const scalar dirY = directionUnit.y();

        if (mag(dirY) > VSMALL)
        {
            const scalar sToMin = (filmYMin_ - entryPoint.y())/dirY;
            const scalar sToMax = (filmYMax_ - entryPoint.y())/dirY;
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
            entryPoint.y() >= filmYMin_
         && entryPoint.y() <= filmYMax_
        )
        {
            filmIntervalStart = 0.0;
            filmIntervalEnd = GREAT;
            haveFilmInterval = true;
        }
    }
    scalar filmEntryOffset = 0.0;

    if (haveFilmInterval && directionMag > VSMALL)
    {
        const scalar dirY = directionUnit.y();

        if (mag(dirY) > VSMALL)
        {
            const scalar targetY = (dirY >= 0) ? filmYMin_ : filmYMax_;
            filmEntryOffset = (targetY - entryPoint.y())/dirY;
            filmEntryOffset = Foam::max(filmEntryOffset, scalar(0));
        }
    }
    forAll(cellCentres, cellI)
    {
        const point& c = cellCentres[cellI];

        if (!searchBox.contains(c))
        {
            continue;
        }

        const bool inFilm = (c.y() >= filmYMin_ && c.y() <= filmYMax_);

        if (!isInBeam(c, axialHalfLength))
        {
            continue;
        }

        if (inFilm) ++metrics.cellsInFilm;
        else        ++metrics.cellsInGas;

        ++metrics.cellsInBeam;
        metrics.totalBeamVolume += cellVolumes[cellI];

        vector r = c - focus_;
        scalar z = (r & direction_);
        r -= z*direction_;
        const scalar R = mag(r);

        const scalar spatialTerm = calculateGaussianIntensity(R);

        const scalar cellAbsorptionCoeff =
            inFilm ? absorptionCoeff_.value() : gasAbsorptionCoeff_.value();

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
            const scalar s =
                directionMag > VSMALL
              ? ((meshPoints[ptI] - entryPoint) & directionUnit)
              : 0.0;
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
        if (directionMag > VSMALL)
        {
            sIn = Foam::max(sIn, scalar(0));
            sOut = Foam::max(sOut, sIn);
        }
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
        scalar transmissionFactor = 1.0;

        if (inFilm)
        {
            const scalar effectiveReflectivity = reflectivity_;
            transmissionFactor = effectiveTransmission(effectiveReflectivity);
        }

        const scalar baseIntensity =
              peakIntensity_.value()
            * temporalAverage
            * spatialTerm
            * transmissionFactor;

        const scalar deltaS = Foam::max(sOut - sIn, VSMALL);

        scalar Ein = baseIntensity;
        scalar Eout = baseIntensity;

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

            if (inFilm) metrics.totalFilmSourceIntegral += cellPower;
            else        metrics.totalGasSourceIntegral  += cellPower;
        }
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
            << "  Max intensity: " << metrics.maxSourceValue/1e12 << " TW/m^3" << nl
            << "  Total power: " << metrics.totalSourceIntegral/1e12 << " TW" << nl
            << "    Film power: " << metrics.totalFilmSourceIntegral/1e12 << " TW" << nl
            << "    Gas power:  " << metrics.totalGasSourceIntegral/1e12 << " TW" << nl
            << "  dt: " << dt << " s" << nl
            << "  Energy this step: " << energyThisStep.value() << " J" << nl
            << "    Film energy this step: " << filmEnergyThisStep << " J" << nl
            << "    Gas energy this step:  " << gasEnergyThisStep  << " J" << nl
            << "  Cumulative energy: " << cumulativeEnergy_ << " J" << nl
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
    cumulativeEnergy_ += stepEnergy;

    if (!continuousLaser_)
    {
        if (depositionActive || trackingPulse_)
        {
            if (!trackingPulse_)
            {
                trackingPulse_ = true;
                pulseEnergyAccumulator_ = 0.0;
                pulseExpectedAccumulator_ = 0.0;
                currentPulseStartTime_ = overlapStart;
            }

            pulseEnergyAccumulator_ += stepEnergy;
            pulseExpectedAccumulator_ += envelope.expectedEnergy;
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
        Info<< "===== LASER DEBUG =====" << nl
            << "Time: " << t << " s (" << t*1e12 << " ps)" << nl
            << "Time step: " << dt << " s (" << dt*1e12 << " ps)" << nl
            << "Laser window: [" << laserStartTime_ << ", " << laserEndTime_ << "] s" << nl
            << "Focus: " << focus_ << nl
            << "  Focus Y: " << focus_.y()*1e6 << " µm" << nl
            << "Film bounds: [" << filmYMin_*1e6 << ", " << filmYMax_*1e6 << "] µm" << nl
            << "Peak intensity: " << peakIntensity_.value() << " W/m²" << nl
            << "Pulse energy: " << pulseEnergy_.value() << " J" << nl
            << "Spot size: " << spotSize_.value()*1e6 << " µm diameter" << nl
            << "Direction: " << direction_ << nl
            << "Absorption coeff: " << absorptionCoeff_.value() << " 1/m" << endl;
    }

    if (timeIndex < lastTimeIndex_)
    {
        cumulativeEnergy_ = 0.0;
        cumulativeFilmEnergy_ = 0.0;
        cumulativeGasEnergy_  = 0.0;
        trackingPulse_ = false;
        pulseEnergyAccumulator_ = 0.0;
        pulseExpectedAccumulator_ = 0.0;
        currentPulseStartTime_ = 0.0;
        pulseCounter_ = 0;
    }
    lastTimeIndex_ = timeIndex;

    const bool singlePulse = (!continuousLaser_ && pulseFrequency_ <= SMALL);

    const scalar effectiveLaserEnd = laserEndTime_;
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
        if (verbose && master && t <= comparisonLaserEnd + 1e-12)
        {
            Info<< "No time overlap - laser inactive this step" << endl;
        }
        finalizePulseEnergyCheck("inactive window", t);
        sourceValid_ = true;
        return;
    }

    const EnvelopeResult envelope =
        evaluateTemporalEnvelope(overlapStart, overlapStart + overlapDuration, dt);

    if (verbose && master && timeIndex % 10 == 0)
    {
        Info<< "Temporal envelope:" << nl
            << "  Active: " << (envelope.active ? "YES" : "NO") << nl
            << "  Temporal integral: " << envelope.temporalIntegral << " s" << nl
            << "  Temporal average: " << envelope.temporalAverage << nl
            << "  Expected energy: " << envelope.expectedEnergy << " J" << endl;
    }

    if (!envelope.active || envelope.temporalAverage <= VSMALL)
    {
        const bool withinActiveWindow =
            (overlapEnd - overlapStart) > SMALL
         && t >= laserStartTime_
         && t <= comparisonLaserEnd;

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

    if (verbose && master && timeIndex % 10 == 0)
    {
        Info<< "Spatial processing results:" << nl
            << "  Cells in beam: " << metrics.cellsInBeam << nl
            << "  Cells in film: " << metrics.cellsInFilm << nl
            << "  Cells in gas: " << metrics.cellsInGas << nl
            << "  Max source value: " << metrics.maxSourceValue << " W/m³" << nl
            << "  Total power: " << metrics.totalSourceIntegral << " W" << nl
            << "  Film power: " << metrics.totalFilmSourceIntegral << " W" << nl
            << "  Limited cells: " << metrics.limitedCells << endl;
    }

    if (metrics.cellsInBeam == 0)
    {
        WarningInFunction
            << "No cells found in laser beam path at time " << t*1e12 << " ps!" << nl
            << "Check focus position, spot size, and mesh resolution." << endl;
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
    return validateParameters() && checkPhysicalBounds();
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
            << "  Max volumetric source: " << maxVolumetricSource_.value()
            << " W/m^3" << nl
            << "  Absorption coefficient: " << absorptionCoeff_.value() << " 1/m" << nl
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
                << "  Cumulative energy: " << cumulativeEnergy_ << " J" << nl
                << "    Film cumulative: " << cumulativeFilmEnergy_ << " J" << nl
                << "    Gas cumulative:  " << cumulativeGasEnergy_  << " J" << endl;
        }
    }
}

} // namespace Foam
