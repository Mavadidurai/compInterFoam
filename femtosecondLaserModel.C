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
        dict.get<scalar>("peakIntensity")
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
    reflectivity_(dict.getOrDefault<scalar>("reflectivity", 0.05)),
    gaussianProfile_(dict.getOrDefault<bool>("gaussianProfile", true)),
    maxReflections_(dict.getOrDefault<label>("maxReflections", 2)),
    continuousLaser_(dict.getOrDefault<bool>("continuousLaser", false)),
    laserStartTime_(dict.getOrDefault<scalar>("laserStartTime", 0.0)),
    laserEndTime_(dict.getOrDefault<scalar>("laserEndTime", 2e-12)),
    sourceValid_(false)
{
    // Normalize direction vector
    direction_ /= mag(direction_) + SMALL;

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
    if (focus_.y() < 8e-6 || focus_.y() > 10e-6) 
    {
        WarningInFunction
            << "Focus Y-coordinate (" << focus_.y()*1e6 << " μm) is outside donor film region (8-10 μm)" << nl
            << "LIFT efficiency may be reduced" << endl;
    }
    
    // Check timing
    scalar pulseDuration = laserEndTime_ - laserStartTime_;
    Info<< "  Pulse duration: " << pulseDuration*1e12 << " ps" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void femtosecondLaserModel::update()
{
    sourceValid_ = false;
}
void femtosecondLaserModel::correct(const scalar time)
{
    // Optional placeholder for time-dependent beam motion/intensity
    // Future implementations can modify focus_, direction_ or
    // peakIntensity_ based on the supplied time value.  For now we simply
    // recompute the source term using the current mesh time information.
    (void)time;  // suppress unused variable warning

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
    const scalar dt = mesh_.time().deltaTValue();
    
    // Simple energy check for pulsed laser
    dimensionedScalar totalEnergy = fvc::domainIntegrate(tSource_()*dt);
    
    // Allow up to 10x pulse energy per timestep (very conservative)
    scalar maxAllowedEnergy = 10.0 * pulseEnergy_.value();
    
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
    const scalar dt = mesh_.time().deltaTValue();
    const label timeIndex = mesh_.time().timeIndex();
    
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
    
    forAll(mesh_.C(), cellI)
    {
        const point& cellCenter = mesh_.C()[cellI];
        
        // Check if in metal film region (8-10μm)
        bool inFilm = (cellCenter.y() >= 8e-6 && cellCenter.y() <= 10e-6);
        if (inFilm) cellsInFilm++;
        
        if (isInBeam(cellCenter))
        {
            cellsInBeam++;
            
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

            if (std::isfinite(intensity) && intensity > 0)
            {
                if (inFilm)
                {
			source[cellI] = intensity;  // No artificial limit!
                }
                else
                {
                    source[cellI] = min(intensity * 0.1, 1e14);  // Reduced outside film
                }
                
                maxSourceValue = max(maxSourceValue, source[cellI]);
                totalSourceIntegral += source[cellI] * mesh_.V()[cellI];
            }
        }
        if (cellsInBeam > 0)
{
    scalar avgIntensityInBeam = totalSourceIntegral / (cellsInBeam * mesh_.V()[0]);
    Info<< "🔍 LASER DIAGNOSTICS:" << nl
        << "  Input peak intensity: " << peakIntensity_.value() << " W/m²" << nl
        << "  Average intensity in beam: " << avgIntensityInBeam << " W/m³" << nl
        << "  Absorption coefficient: " << absorptionCoeff_.value() << " 1/m" << nl
        << "  Spot radius: " << spotSize_.value()/2.0*1e6 << " μm" << nl
        << "  Beam area: " << 3.14159*sqr(spotSize_.value()/2.0)*1e12 << " μm²" << endl;
}
    }

    // Reduce across processors
    reduce(cellsInBeam, sumOp<label>());
    reduce(cellsInFilm, sumOp<label>());
    reduce(maxSourceValue, maxOp<scalar>());
    reduce(totalSourceIntegral, sumOp<scalar>());

    sourceValid_ = true;

    // Report laser activity
    if (laserActive && (timeIndex % 10 == 0))
    {
        dimensionedScalar totalEnergyDeposited = 
            fvc::domainIntegrate(source * dt);
            
        Info<< "🔥 LASER ENERGY DEPOSITION:" << nl
            << "  Time: " << t*1e12 << " ps" << nl
            << "  Temporal factor: " << temporalTerm << nl
            << "  Cells in beam: " << cellsInBeam << nl
            << "  Cells in metal film: " << cellsInFilm << nl
            << "  Max intensity: " << maxSourceValue/1e12 << " TW/m³" << nl
            << "  Total power: " << totalSourceIntegral/1e12 << " TW" << nl
            << "  Energy this step: " << totalEnergyDeposited.value()/1e-6 << " μJ" << endl;
            
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
}

} // End namespace Foam

// ************************************************************************* //
