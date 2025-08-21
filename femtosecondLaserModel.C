/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom LIFT Solver
    
    Description
    Femtosecond laser model for LIFT (Laser-Induced Forward Transfer) simulation.
    Handles substrate-backed geometry with realistic laser-material interaction.

\*---------------------------------------------------------------------------*/

#include "femtosecondLaserModel.H"
#include "fvc.H"
#include "fvm.H"
#include "mathematicalConstants.H"
#include "HashSet.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(femtosecondLaserModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

femtosecondLaserModel::femtosecondLaserModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    peakIntensity_
    (
        "peakIntensity",
        dimPower/dimArea,
        dict.getOrDefault<scalar>("peakIntensity", 0.0)
    ),
    pulseWidth_
    (
        "pulseWidth",
        dimTime,
        dict.get<scalar>("pulseWidth")
    ),
    wavelength_
    (
        "wavelength",
        dimLength,
        dict.get<scalar>("wavelength")
    ),
    absorptionCoeff_
    (
        "absorptionCoeff",
        dimless/dimLength,
        dict.getOrDefault<scalar>("absorptionCoeff", 5e6)
    ),
    spotSize_
    (
        "spotSize",
        dimLength,
        dict.get<scalar>("spotSize")
    ),
    pulseEnergy_
    (
        "pulseEnergy",
        dimEnergy,
        dict.get<scalar>("pulseEnergy")
    ),
    direction_(dict.get<vector>("direction")),
    focus_(dict.get<point>("focus")),
    initialFocus_(focus_),
    scanVelocity_(dict.getOrDefault<vector>("scanVelocity", vector::zero)),
    pulseFrequency_(dict.getOrDefault<scalar>("pulseFrequency", 0.0)),
    pulseDutyCycle_(dict.getOrDefault<scalar>("pulseDutyCycle", 1.0)),
    reflectivity_(dict.getOrDefault<scalar>("reflectivity", 0.05)),
    gaussianProfile_(dict.getOrDefault<bool>("gaussianProfile", true)),
    maxReflections_(dict.getOrDefault<label>("maxReflections", 2)),
    continuousLaser_(dict.getOrDefault<bool>("continuousLaser", false)),
    laserStartTime_(dict.getOrDefault<scalar>("laserStartTime", 0.0)),
    laserEndTime_(dict.getOrDefault<scalar>("laserEndTime", 2e-12)),
    filmYMin_(dict.getOrDefault<scalar>("filmYMin", 8e-6)),
    filmYMax_(dict.getOrDefault<scalar>("filmYMax", 10e-6)),
    sourceValid_(false),
    cumulativeEnergy_(0.0),
    lastTimeIndex_(mesh.time().timeIndex())
{
    // Normalize direction vector
    direction_ /= mag(direction_) + SMALL;

    // Derive peak intensity from pulse energy if not provided
    scalar derivedPeak =
        pulseEnergy_.value()
        /
        (
            constant::mathematical::pi*sqr(spotSize_.value()/2.0)*
            pulseWidth_.value()
        );

    if (dict.found("peakIntensity"))
    {
        scalar diff = mag(peakIntensity_.value() - derivedPeak);
        scalar tol = 0.05*derivedPeak;
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
            "peakIntensity",
            dimPower/dimArea,
            derivedPeak
        );
    }

    // Validate parameters
    if (!valid())
    {
        FatalErrorInFunction
            << "Invalid laser parameters" << nl
            << "  Peak intensity: " << peakIntensity_.value() << nl
            << "  Pulse width: " << pulseWidth_.value() << nl
            << "  Wavelength: " << wavelength_.value() << nl
            << "  Spot size: " << spotSize_.value()
            << abort(FatalError);
    }

    Info<< "✅ Femtosecond laser model initialized:" << nl
        << "  Mode: " << (continuousLaser_ ? "Continuous" : "Pulsed") << nl
        << "  Peak intensity: " << peakIntensity_.value() << " W/m²" << nl
        << "  Pulse width: " << pulseWidth_.value() << " s" << nl
        << "  Wavelength: " << wavelength_.value() << " m" << nl
        << "  Spot size: " << spotSize_.value() << " m" << nl
        << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
        << "  Focus: " << focus_ << nl
        << "  Direction: " << direction_ << nl
        << "  Active time: " << laserStartTime_ << " to " << laserEndTime_ << " s" << endl;
        
    // Validate focus position for LIFT geometry
    if (focus_.y() < filmYMin_ || focus_.y() > filmYMax_)
    {
        WarningInFunction
            << "Focus Y-coordinate (" << focus_.y()*1e6
            << " μm) is outside donor film region ("
            << filmYMin_*1e6 << "-" << filmYMax_*1e6 << " μm)" << nl
            << "LIFT efficiency may be reduced" << endl;
    }
    
    // Check timing
    scalar pulseDuration = laserEndTime_ - laserStartTime_;
    Info<< "  Pulse duration: " << pulseDuration*1e12 << " ps" << endl;
    
    // Warn about any unhandled dictionary entries
    wordHashSet handled
    (
        wordList
        {
            "peakIntensity",
            "pulseWidth",
            "wavelength",
            "absorptionCoeff",
            "spotSize",
            "pulseEnergy",
            "direction",
            "focus",
            "scanVelocity",
            "pulseFrequency",
            "pulseDutyCycle",
            "reflectivity",
            "gaussianProfile",
            "maxReflections",
            "continuousLaser",
            "laserStartTime",
            "laserEndTime",
            "filmYMin",
            "filmYMax"
        }
    );

    forAllConstIter(dictionary, dict_, iter)
    {
        const word& key = iter().keyword();
        if (!handled.found(key))
        {
            WarningInFunction
                << "Ignoring unhandled entry '" << key << "'" << endl;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void femtosecondLaserModel::update()
{
        // Handle laser scanning by updating the focal point based on
    // the specified scan velocity.  The initial focus position is
    // stored separately so that the motion is always measured from
    // the original location.
    const Foam::scalar currentTime = mesh_.time().value();

    
    if (mag(scanVelocity_) > VSMALL)
    {
focus_ = initialFocus_ + scanVelocity_ * mesh_.time().value();
    }

    // Determine if the laser pulse is active at the supplied time.  For a
    // continuous laser this always evaluates to true inside the time window
    // [laserStartTime_, laserEndTime_].  For pulsed operation the duty cycle
    // and repetition frequency are honoured.
    bool pulseActive = true;

    if (!continuousLaser_)
    {
        pulseActive = (currentTime >= laserStartTime_ && currentTime <= laserEndTime_);

        if (pulseActive && pulseFrequency_ > SMALL)
        {
            const scalar period = 1.0/pulseFrequency_;
            const scalar localTime = currentTime - laserStartTime_;
            const scalar timeInPeriod = std::fmod(localTime, period);
            pulseActive = (timeInPeriod <= pulseDutyCycle_*period);
        }
    }

    if (!pulseActive)
    {
        // Ensure the source field exists and set it to zero so downstream
        // code can safely access it without additional checks.
        if (!tSource_.valid())
        {
            tSource_ = tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "laserSource",
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
        else
        {
            tSource_.ref() = dimensionedScalar("zero", dimPower/dimVolume, 0.0);
        }

        sourceValid_ = true;
        return;
    }
    sourceValid_ = false;
}
void femtosecondLaserModel::correct(const scalar currentTime)
{
    (void)currentTime;  // suppress unused variable warning
    sourceValid_ = false;
    calculateSource();
}

bool femtosecondLaserModel::validateParameters() const
{
    bool valid = true;

    // Check dimensions
    if (peakIntensity_.dimensions() != dimPower/dimArea ||
        pulseWidth_.dimensions() != dimTime ||
        wavelength_.dimensions() != dimLength ||
        absorptionCoeff_.dimensions() != dimless/dimLength ||
        spotSize_.dimensions() != dimLength ||
        pulseEnergy_.dimensions() != dimEnergy)
    {
        valid = false;
    }

    // Check values
    if (peakIntensity_.value() <= 0 ||
        pulseWidth_.value() <= 0 ||
        wavelength_.value() <= 0 ||
        spotSize_.value() <= 0 ||
        pulseEnergy_.value() <= 0)
    {
        valid = false;
    }

    return valid;
}

bool femtosecondLaserModel::checkPhysicalBounds() const
{
    bool valid = true;

    // Check femtosecond regime (10fs - 10ps acceptable for debug)
    if (!continuousLaser_ && (pulseWidth_.value() < 1e-16 || pulseWidth_.value() > 1e-11))
    {
        WarningInFunction
            << "Pulse width outside reasonable range: "
            << pulseWidth_.value() << " s" << endl;
    }

    // Check wavelength (visible to near-IR)
    if (wavelength_.value() < 1e-7 || wavelength_.value() > 2e-6)
    {
        WarningInFunction
            << "Wavelength outside typical range: "
            << wavelength_.value() << " m" << endl;
    }

    // Check intensity (reasonable for LIFT)
    if (peakIntensity_.value() > 1e16)
    {
        WarningInFunction
            << "Very high peak intensity: "
            << peakIntensity_.value() << " W/m²" << endl;
    }

    return valid;
}

scalar femtosecondLaserModel::calculateGaussianIntensity
(
    const scalar R,
    const scalar z
) const
{
    if (gaussianProfile_)
    {
        // Gaussian beam profile: I = I0 * exp(-2*r²/w²)
        scalar beamRadius = spotSize_.value() / 2.0;  // Convert diameter to radius
        return exp(-2.0*sqr(R)/sqr(beamRadius));
    }
    else
    {
        // Top-hat profile
        return R <= (spotSize_.value() / 2.0) ? 1.0 : 0.0;
    }
}

bool femtosecondLaserModel::isInBeam(const point& p) const
{
    // Calculate distance from focus
    vector r = p - focus_;
    scalar z = (r & direction_);
    r -= z*direction_;
    scalar R = mag(r);
    
    // Check if within beam radius (with some tolerance)
    scalar beamRadius = spotSize_.value() / 2.0;
    scalar maxRadius = beamRadius * 3.0;  // 3-sigma cutoff for Gaussian
    
    // Check z-direction limits (allow some propagation)
    scalar maxZ = max(spotSize_.value(), 2e-6);  // At least 2 μm
    
    return (R <= maxRadius) && (z >= -maxZ) && (z <= maxZ);
}

bool femtosecondLaserModel::checkEnergyConservation() const
{
    if (!tSource_.valid() || continuousLaser_)
    {
        return true;
    }

    const scalar t = mesh_.time().value();
    (void)t;
    const dimensionedScalar dt = mesh_.time().deltaT();
    
       // Simple energy check for pulsed laser
    dimensionedScalar totalEnergy = fvc::domainIntegrate(tSource_()*dt);

    // Allow up to 10x energy implied by peak intensity per timestep (very conservative)
    dimensionedScalar expectedEnergy
    (
        "expectedEnergy",
        dimEnergy,
        peakIntensity_.value()
       *constant::mathematical::pi*sqr(spotSize_.value()/2.0)
       *pulseWidth_.value()
    );

    scalar maxAllowedEnergy = 10.0 * expectedEnergy.value();
    
    return totalEnergy.value() <= maxAllowedEnergy;
}

void femtosecondLaserModel::calculateSource() const
{
    if (sourceValid_)
    {
        return;
    }

    // Create source field if needed
    if (!tSource_.valid())
    {
        tSource_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "laserSource",
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
    
    // Get current time
    const scalar t = mesh_.time().value();
    const dimensionedScalar dt = mesh_.time().deltaT();
    const label timeIndex = mesh_.time().timeIndex();
    
    // Reset cumulative energy if simulation restarted
    if (timeIndex < lastTimeIndex_)
    {
        cumulativeEnergy_ = 0.0;
    }
    lastTimeIndex_ = timeIndex;
    
    // Check if laser is active
    bool laserActive = false;
    
    if (continuousLaser_)
    {
        laserActive = (t >= laserStartTime_ && t <= laserEndTime_);
    }
    else
    {
        // For pulsed laser, check if within pulse duration
        laserActive = (t >= laserStartTime_ && t <= laserEndTime_);
    }
    
    // Report laser status
    if (timeIndex % 50 == 0)  // Every 50 timesteps
    {
        Info<< "🔍 LASER @ step=" << timeIndex 
            << ": t=" << t*1e12 << "ps, " 
            << (laserActive ? "ACTIVE" : "inactive") 
            << " (window: " << laserStartTime_*1e12 << "-" << laserEndTime_*1e12 << "ps)" << endl;
    }

    if (!laserActive)
    {
        sourceValid_ = true;
        return;
    }

    // Calculate temporal profile
    scalar temporalTerm = 1.0;
    if (!continuousLaser_)
    {
        // Gaussian temporal profile centered in the pulse
        scalar pulseDuration = laserEndTime_ - laserStartTime_;
        scalar pulseCenter = laserStartTime_ + pulseDuration/2.0;
        scalar timeFromCenter = t - pulseCenter;
        scalar sigma = pulseWidth_.value() / (2.0 * sqrt(2.0 * log(2.0)));  // FWHM to sigma
        
        temporalTerm = exp(-0.5 * sqr(timeFromCenter/sigma));
    }

    // Apply laser heating
    label cellsInBeam = 0;
    label cellsInFilm = 0;
    scalar maxSourceValue = 0.0;
    scalar totalSourceIntegral = 0.0;
    scalar totalBeamVolume = 0.0;

    
    forAll(mesh_.C(), cellI)
    {
        const point& cellCenter = mesh_.C()[cellI];
        
        // Check if in metal film region
        bool inFilm = (cellCenter.y() >= filmYMin_ && cellCenter.y() <= filmYMax_);
        if (inFilm) cellsInFilm++;
        
        if (isInBeam(cellCenter))
        {
            cellsInBeam++;
            totalBeamVolume += mesh_.V()[cellI];

            
            vector r = cellCenter - focus_;
            scalar z = (r & direction_);
            r -= z*direction_;
            scalar R = mag(r);
            
            scalar spatialTerm = calculateGaussianIntensity(R, z);
            scalar absorptionTerm = exp(-absorptionCoeff_.value() * max(z, 0.0));

            scalar intensity = peakIntensity_.value() *
                              temporalTerm *
                              spatialTerm *
                              absorptionTerm *
                              (1.0 - reflectivity_);

            scalar volumetricIntensity = intensity * absorptionCoeff_.value();

            if (std::isfinite(volumetricIntensity) && volumetricIntensity > 0)
            {
                if (inFilm)
                {
                    source[cellI] = volumetricIntensity;  // No artificial limit!
                }
                else
                {
                    source[cellI] = min(volumetricIntensity * 0.1, 1e14);  // Reduced outside film
                }
                
                maxSourceValue = max(maxSourceValue, source[cellI]);
                totalSourceIntegral += source[cellI] * mesh_.V()[cellI];
            }
        }

    }

    // Reduce across processors
    reduce(cellsInBeam, sumOp<label>());
    reduce(cellsInFilm, sumOp<label>());
    reduce(maxSourceValue, maxOp<scalar>());
    reduce(totalSourceIntegral, sumOp<scalar>());
    reduce(totalBeamVolume, sumOp<scalar>());

    scalar avgIntensityInBeam = 0.0;
    if (totalBeamVolume > 0)
    {
        avgIntensityInBeam = totalSourceIntegral / totalBeamVolume;
    }
    dimensionedScalar totalEnergyDeposited =
        fvc::domainIntegrate(source * dt);
    cumulativeEnergy_ += totalEnergyDeposited.value();

    if (!checkEnergyConservation())
    {
        WarningInFunction
            << "Energy conservation check failed at timestep "
            << timeIndex << ": energy = "
            << totalEnergyDeposited.value() << " J" << endl;
    }

    Info<< "🔍 LASER DIAGNOSTICS:" << nl
        << "  Input peak intensity: " << peakIntensity_.value() << " W/m²" << nl
        << "  Average volumetric intensity in beam: " << avgIntensityInBeam << " W/m³" << nl
        << "  Absorption coefficient: " << absorptionCoeff_.value() << " 1/m" << nl
        << "  Spot radius: " << spotSize_.value()/2.0*1e6 << " μm" << nl
        << "  Beam area: " << 3.14159*sqr(spotSize_.value()/2.0)*1e12 << " μm²" << endl;

    sourceValid_ = true;

    // Report laser activity
    if (laserActive && (timeIndex % 10 == 0))
    {
        Info<< "🔥 LASER ENERGY DEPOSITION:" << nl
            << "  Time: " << t*1e12 << " ps" << nl
            << "  Temporal factor: " << temporalTerm << nl
            << "  Cells in beam: " << cellsInBeam << nl
            << "  Cells in metal film: " << cellsInFilm << nl
            << "  Max intensity: " << maxSourceValue/1e12 << " TW/m³" << nl
            << "  Total power: " << totalSourceIntegral/1e12 << " TW" << nl
            << "  Time step: " << dt.value() << " s" << nl
            << "  Energy this step: " << totalEnergyDeposited.value() << " J" << nl
            << "  Cumulative energy: " << cumulativeEnergy_ << " J" << endl;

        // Diagnostics
        if (cellsInBeam == 0)
        {
            WarningInFunction
                << "No cells in laser beam!" << nl
                << "Focus: " << focus_ << nl
                << "Domain bounds: " << mesh_.bounds() << endl;
        }

        if (cellsInFilm == 0 && cellsInBeam > 0)
        {
            WarningInFunction
                << "Beam hits " << cellsInBeam
                << " cells but none in metal film region" << endl;
        }
    }
}

tmp<volScalarField> femtosecondLaserModel::source() const
{
    calculateSource();
    return tSource_;
}

bool femtosecondLaserModel::valid() const
{
    return validateParameters() && checkPhysicalBounds();
}

void femtosecondLaserModel::write() const
{
    Info<< "Femtosecond laser model status:" << nl
        << "  Mode: " << (continuousLaser_ ? "Continuous" : "Pulsed") << nl
        << "  Peak intensity: " << peakIntensity_.value() << " W/m²" << nl
        << "  Pulse width: " << pulseWidth_.value() << " s" << nl
        << "  Wavelength: " << wavelength_.value() << " m" << nl
        << "  Spot size: " << spotSize_.value() << " m" << nl
        << "  Pulse energy: " << pulseEnergy_.value() << " J" << nl
        << "  Absorption coefficient: " << absorptionCoeff_.value() << " 1/m" << nl
        << "  Reflectivity: " << reflectivity_ << nl
        << "  Focus: " << focus_ << nl
        << "  Direction: " << direction_ << nl
        << "  Active time: " << laserStartTime_ << " to " << laserEndTime_ << " s" << endl;

    if (tSource_.valid())
    {
        const dimensionedScalar maxIntensity = max(tSource_());
        const dimensionedScalar totalEnergy =
            fvc::domainIntegrate(tSource_() * mesh_.time().deltaT());

        Info<< "Source statistics:" << nl
            << "  Maximum intensity: " << maxIntensity.value() << " W/m³" << nl
            << "  Energy this timestep: " << totalEnergy.value() << " J" << endl;
    }

    Info<< "  Cumulative energy: " << cumulativeEnergy_ << " J" << endl;
}
} // End namespace Foam

// ************************************************************************* //
