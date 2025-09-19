/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Femtosecond laser model for LIFT (implementation)
\*---------------------------------------------------------------------------*/
#include "femtosecondLaserModel.H"
#include "fvc.H"
#include "fvm.H"
#include "mathematicalConstants.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "treeBoundBox.H"
#include "treeDataCell.H"

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
            dict.getOrDefault<scalar>("maxVolumetricSource", 1e18)
        )
    ),    
    direction_(dict.get<vector>("direction")),
    focus_(dict.get<point>("focus")),
    initialFocus_(focus_),
    scanVelocity_(dict.getOrDefault<vector>("scanVelocity", vector::zero)),
    pulseFrequency_(dict.getOrDefault<scalar>("pulseFrequency", 0.0)),
    pulseDutyCycle_(dict.getOrDefault<scalar>("pulseDutyCycle", 1.0)),
    reflectivity_(dict.getOrDefault<scalar>("reflectivity", 0.05)),
    transmission_(dict.getOrDefault<scalar>("transmission", -1.0)),
    incidenceAngle_(dict.getOrDefault<scalar>("incidenceAngle", 0.0)),
    gaussianProfile_(dict.getOrDefault<bool>("gaussianProfile", true)),
    continuousLaser_(dict.getOrDefault<bool>("continuousLaser", false)),
    laserStartTime_(dict.getOrDefault<scalar>("laserStartTime", 0.0)),
    laserEndTime_(dict.getOrDefault<scalar>("laserEndTime", 2e-12)),
    filmYMin_(dict.getOrDefault<scalar>("filmYMin", 8e-6)),
    filmYMax_(dict.getOrDefault<scalar>("filmYMax", 10e-6)),
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
    trackingPulse_(false)
{
    // normalize direction
    direction_ /= mag(direction_) + SMALL;
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
            WarningInFunction
                << "Supplied peakIntensity (" << peakIntensity_.value()
                << " W/m^2) differs from value derived from pulseEnergy "
                << derivedPeak << " W/m^2" << endl;
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

    if (verbose)
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
    const scalar expected =
        dict_.getOrDefault<scalar>("filmThicknessExpected", 7.14e-8);
    const scalar tol =
        dict_.getOrDefault<scalar>("filmThicknessTolerance", 0.1*expected);
    if (filmThickness <= 0)
    {
        WarningInFunction
            << "Non-positive film thickness: " << filmThickness
            << " (check filmYMin, filmYMax)" << endl;
    }
    if (mag(filmThickness - expected) > tol)
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
    if (peakIntensity_.value() > 6e16)
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
bool femtosecondLaserModel::isInBeam(const point& p) const
{
    vector r = p - focus_;
    scalar z = (r & direction_);
    r -= z*direction_;
    const scalar R = mag(r);

    const scalar beamRadius = spotSize_.value()/2.0;
    const scalar maxRadius  = 3.0*beamRadius;         // ~3-sigma
    const scalar maxZ       = max(spotSize_.value(), 2e-6);

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
    const scalar configured = pulseEnergy_.value();
    const scalar reference = max(max(expected, configured), scalar(1e-12));
    const scalar tolerance = max(scalar(1e-12), 0.01*reference);
    const scalar diffExpected = mag(pulseEnergyAccumulator_ - expected);
    const scalar diffConfigured = mag(pulseEnergyAccumulator_ - configured);
    const scalar expectedMismatch = mag(expected - configured);

    ++pulseCounter_;

    if (verbose)
    {
        Info<< "LASER PULSE ENERGY CHECK:" << nl
            << "  Pulse index:      " << pulseCounter_ << nl
            << "  Context:          " << context << nl
            << "  Start time:       " << currentPulseStartTime_ << " s" << nl
            << "  End time:         " << currentTime << " s" << nl
            << "  Deposited energy:      " << pulseEnergyAccumulator_ << " J" << nl
            << "  Expected (integrated): " << expected << " J" << nl
            << "  Configured pulse:     " << configured << " J" << nl
            << "  Diff vs expected:     " << diffExpected << " J" << nl
            << "  Diff vs configured:   " << diffConfigured << " J" << nl
            << "  Expected-config diff: " << expectedMismatch << " J" << nl
            << "  Tolerance:            " << tolerance << " J" << endl;
    }

    if (diff > tolerance)
    {
        WarningInFunction
            << "Laser pulse energy mismatch after pulse " << pulseCounter_
            << ": deposited " << pulseEnergyAccumulator_ << " J vs requested "
            << configured << " J (tolerance " << tolerance << ")" << endl;
    }
    else if (expectedMismatch > tolerance && verbose)
    {
        Info<< "  Note: integrated expectation differs from configured energy by "
            << expectedMismatch << " J" << endl;
    }

    trackingPulse_ = false;
    pulseEnergyAccumulator_ = 0.0;
    pulseExpectedAccumulator_ = 0.0;    
    currentPulseStartTime_ = 0.0;
}
//------------------------------------------------------------------------------
// source compute (with pulsed window + duty cycle + temporal Gaussian)
//------------------------------------------------------------------------------
void femtosecondLaserModel::calculateSource() const
{
    if (sourceValid_) return;

    // Allocate/reset field
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

    // Time bookkeeping
    const scalar t = mesh_.time().value();
    const dimensionedScalar dtDim = mesh_.time().deltaT();
    const scalar dt = dtDim.value();
    const scalar tStart = t - dt;
    const label timeIndex = mesh_.time().timeIndex();
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

    const scalar overlapStart = max(tStart, laserStartTime_);
    const scalar overlapEnd   = min(t, laserEndTime_);
    const scalar overlapDuration = overlapEnd - overlapStart;

    if (overlapDuration <= VSMALL)
    {
        finalizePulseEnergyCheck("inactive window", t);
        sourceValid_ = true; // keep zero field
        return;
    }

    scalar temporalIntegral = 0.0;
    scalar expectedEnergyThisStep = 0.0;

    const auto gaussianIntegral =
        [](const scalar a, const scalar b, const scalar center, const scalar sigma)
        {
            const scalar invSqrt2Sigma = 1.0/(sqrt(2.0)*sigma);
            const scalar prefactor = sigma*sqrt(constant::mathematical::pi/2.0);
            return prefactor
                * (std::erf((b - center)*invSqrt2Sigma)
                 - std::erf((a - center)*invSqrt2Sigma));
        };

    if (continuousLaser_)
    {
        temporalIntegral = overlapDuration;
    }
    else if (pulseFrequency_ > SMALL)
    {
        // Repeated fs pulses with a Gaussian inside each ON-window
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
            const scalar integral = gaussianIntegral(clipStart, clipEnd, center, sigma);

            if (integral > VSMALL)
            {
                temporalIntegral += integral;
            }
        }

        if (temporalIntegral > VSMALL)
        {
            expectedEnergyThisStep =
                pulseEnergy_.value()*temporalIntegral/max(fullPulseIntegral, VSMALL);
        }
    }
    else
    {
        // Single Gaussian across the entire active window
        const scalar sigma  = pulseWidth_.value()/(2.0*sqrt(2.0*log(2.0)));
        const scalar center = 0.5*(laserStartTime_ + laserEndTime_);

        temporalIntegral = gaussianIntegral(overlapStart, overlapEnd, center, sigma);

        if (temporalIntegral > VSMALL)
        {
            const scalar fullIntegral = sigma*sqrt(2.0*constant::mathematical::pi);
            expectedEnergyThisStep =
                pulseEnergy_.value()*temporalIntegral/max(fullIntegral, VSMALL);
        }
    }

    const bool laserActive = temporalIntegral > VSMALL;

    scalar temporalAverage = 0.0;
    if (laserActive)
    {
        temporalAverage = temporalIntegral/max(dt, VSMALL);
        temporalAverage = min(scalar(1.0), max(temporalAverage, scalar(0.0)));
    }

    // Negligible envelope → skip work this step
    if (!laserActive || temporalAverage <= VSMALL)
    {
        finalizePulseEnergyCheck("pulse window complete", t);
        sourceValid_ = true;
        return;
    }
    // --- Apply laser heating over cells ---
    const scalar beamRadius = spotSize_.value()/2.0;
    const scalar radialHalfWidth = 3.0*beamRadius; // ~3-sigma laterally

    scalar axialHalfLength = max(spotSize_.value(), 2e-6);
    if (absorptionCoeff_.value() > VSMALL)
    {
        const scalar absorptionDepth = 1.0/absorptionCoeff_.value();
        axialHalfLength = max(axialHalfLength, 3.0*absorptionDepth);
    }

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

    labelList treeCandidates = mesh_.cellTree().findBox(searchBox);

    DynamicList<label> candidateCells;
    candidateCells.reserve(treeCandidates.size());

    forAll(treeCandidates, idx)
    {
        const label cellI = treeCandidates[idx];
        const point& c = mesh_.C()[cellI];

        if (searchBox.contains(c))
        {
            candidateCells.append(cellI);
        }
    }

    if (candidateCells.empty())
    {
        forAll(mesh_.C(), cellI)
        {
            const point& c = mesh_.C()[cellI];

            if (searchBox.contains(c))
            {
                candidateCells.append(cellI);
            }
        }
    }

    candidateCells.shrink();

    label cellsInBeam = 0, cellsInFilm = 0, cellsInGas = 0;
    scalar maxSourceValue = 0.0;
    scalar totalSourceIntegral = 0.0;
    scalar totalBeamVolume = 0.0;
    scalar totalFilmSourceIntegral = 0.0;
    scalar totalGasSourceIntegral  = 0.0;
    label limitedCells = 0;

    const scalar maxSourceCap = maxVolumetricSource_.value();
    const bool limitLaserSource = maxSourceCap > SMALL;
    if (transmission_ >= 0)
    {
        WarningInFunction
            << "transmission overrides reflectivity" << endl;
    }

    forAll(candidateCells, candidateI)
    {
        const label cellI = candidateCells[candidateI];
        const point& c = mesh_.C()[cellI];

        const bool inFilm = (c.y() >= filmYMin_ && c.y() <= filmYMax_);
        if (!isInBeam(c)) continue;
        if (inFilm) ++cellsInFilm;
        else        ++cellsInGas;

        ++cellsInBeam;
        totalBeamVolume += mesh_.V()[cellI];

        vector r = c - focus_;
        scalar z = (r & direction_);
        r -= z*direction_;
        const scalar R = mag(r);

        const scalar spatialTerm = calculateGaussianIntensity(R);

        const scalar cellAbsorptionCoeff =
            inFilm ? absorptionCoeff_.value() : gasAbsorptionCoeff_.value();

        const scalar absorptionTerm =
            (cellAbsorptionCoeff > VSMALL)
          ? exp(-cellAbsorptionCoeff*max(z, 0.0))
          : 1.0;

        scalar transmissionFactor = 1.0 - reflectivity_;
        if (transmission_ >= 0) transmissionFactor = transmission_;
        else if (incidenceAngle_ > VSMALL)
        {
            const scalar n1 = 1.0;
            const scalar sqrtR = sqrt(max(0.0, reflectivity_));
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

        const scalar intensity =
              peakIntensity_.value()
            * temporalAverage
            * spatialTerm
            * absorptionTerm
            * transmissionFactor;

        const scalar volIntensity = intensity * cellAbsorptionCoeff;

        if (std::isfinite(volIntensity) && volIntensity > 0)
        {
            scalar limitedValue = volIntensity;

            if (limitLaserSource && limitedValue > maxSourceCap)
            {
                limitedValue = maxSourceCap;
                ++limitedCells;
            }

            limitedValue = Foam::max(limitedValue, scalar(0));

            source[cellI] = limitedValue;


            const scalar cellPower = source[cellI] * mesh_.V()[cellI];

            maxSourceValue       = max(maxSourceValue, source[cellI]);
            totalSourceIntegral += cellPower;

            if (inFilm) totalFilmSourceIntegral += cellPower;
            else        totalGasSourceIntegral  += cellPower;
        }
    }

    // parallel reductions
    reduce(cellsInBeam, sumOp<label>());
    reduce(cellsInFilm, sumOp<label>());
    reduce(cellsInGas,  sumOp<label>());
    reduce(maxSourceValue, maxOp<scalar>());
    reduce(totalSourceIntegral, sumOp<scalar>());
    reduce(totalBeamVolume, sumOp<scalar>());
    reduce(totalFilmSourceIntegral, sumOp<scalar>());
    reduce(totalGasSourceIntegral, sumOp<scalar>());
    reduce(limitedCells, sumOp<label>());

    if (limitLaserSource && limitedCells > 0 && verbose)
    {
        Info<< "Laser source limited to " << maxSourceCap
            << " W/m^3 in " << limitedCells << " cells" << endl;
    }
    const scalar avgIntensityInBeam =
        (totalBeamVolume > 0) ? totalSourceIntegral/totalBeamVolume : 0.0;

    const dimensionedScalar energyThisStep =
        fvc::domainIntegrate(source * dtDim);
    const scalar stepEnergy = energyThisStep.value();
    cumulativeEnergy_ += stepEnergy;

    if (!continuousLaser_)
    {
        if (laserActive || trackingPulse_)
        {
            if (!trackingPulse_)
            {
                trackingPulse_ = true;
                pulseEnergyAccumulator_ = 0.0;
                pulseExpectedAccumulator_ = 0.0;
                currentPulseStartTime_ = overlapStart;
            }

            pulseEnergyAccumulator_ += stepEnergy;
            pulseExpectedAccumulator_ += expectedEnergyThisStep;
        }
    }
    const scalar filmEnergyThisStep = totalFilmSourceIntegral * dt;
    const scalar gasEnergyThisStep  = totalGasSourceIntegral  * dt;

    cumulativeFilmEnergy_ += filmEnergyThisStep;
    cumulativeGasEnergy_  += gasEnergyThisStep;

    if (!checkEnergyConservation())
    {
        WarningInFunction
            << "Energy check failed at step " << timeIndex
            << " (E=" << energyThisStep.value() << " J)" << endl;
    }

    if (verbose)
    {
        Info<< "LASER DIAGNOSTICS:" << nl
            << "  Input peak intensity: " << peakIntensity_.value() << " W/m^2" << nl
            << "  Average volumetric intensity in beam: " << avgIntensityInBeam << " W/m^3" << nl
            << "  Absorption coefficient (film): " << absorptionCoeff_.value() << " 1/m" << nl
            << "  Absorption coefficient (gas):  " << gasAbsorptionCoeff_.value() << " 1/m" << nl
            << "  Spot radius: " << spotSize_.value()/2.0*1e6 << " µm" << nl
            << "  Beam area: "
            << constant::mathematical::pi*sqr(spotSize_.value()/2.0)*1e12
            << " µm^2" << endl;

        if (laserActive && (timeIndex % 10 == 0))
        {
            Info<< "LASER ENERGY DEPOSITION:" << nl
                << "  Time: " << t*1e12 << " ps" << nl
                << "  Temporal factor: " << temporalAverage << nl
                << "  Cells in beam: " << cellsInBeam << nl
                << "  Cells in metal film: " << cellsInFilm << nl
                << "  Cells in gas: " << cellsInGas << nl
                << "  Max intensity: " << maxSourceValue/1e12 << " TW/m^3" << nl
                << "  Total power: " << totalSourceIntegral/1e12 << " TW" << nl
                << "    Film power: " << totalFilmSourceIntegral/1e12 << " TW" << nl
                << "    Gas power:  " << totalGasSourceIntegral/1e12 << " TW" << nl
                << "  dt: " << dt << " s" << nl
                << "  Energy this step: " << energyThisStep.value() << " J" << nl
                << "    Film energy this step: " << filmEnergyThisStep << " J" << nl
                << "    Gas energy this step:  " << gasEnergyThisStep  << " J" << nl
                << "  Cumulative energy: " << cumulativeEnergy_ << " J" << nl
                << "    Cumulative film: " << cumulativeFilmEnergy_ << " J" << nl
                << "    Cumulative gas:  " << cumulativeGasEnergy_  << " J" << endl;

            if (cellsInBeam == 0)
            {
                WarningInFunction
                    << "No cells in laser beam!"
                    << " Focus: " << focus_
                    << " Bounds: " << mesh_.bounds() << endl;
            }
            if (cellsInFilm == 0 && cellsInBeam > 0)
            {
                WarningInFunction
                    << "Beam hits " << cellsInBeam
                    << " cells but none in metal film region" << endl;
            }
            else if (cellsInGas > 0 && gasAbsorptionCoeff_.value() <= VSMALL)
            {
                Info<< "    Gas absorption disabled (gasAbsorptionCoeff ≈ 0);"
                    << " no direct heating deposited off-film." << endl;
            }            
        }
    }

    sourceValid_ = true;
}

//------------------------------------------------------------------------------
tmp<volScalarField> femtosecondLaserModel::source() const
{
    calculateSource();
    return tSource_;
}

//------------------------------------------------------------------------------
bool femtosecondLaserModel::valid() const
{
    return validateParameters() && checkPhysicalBounds();
}

//------------------------------------------------------------------------------
void femtosecondLaserModel::write() const
{
    if (verbose)
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

