/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024

    Header: femtosecond laser model for LIFT
\*---------------------------------------------------------------------------*/
#ifndef femtosecondLaserModel_H
#define femtosecondLaserModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "mathematicalConstants.H"
#include "tmp.H"

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
    dimensionedScalar   spotSize_;          // [m] (diameter)
    dimensionedScalar   pulseEnergy_;       // [J]

    // Geometry / timing
    vector              direction_;         // unit vector
    point               focus_;
    point               initialFocus_;
    vector              scanVelocity_;      // [m/s]
    scalar              pulseFrequency_;    // [1/s]
    scalar              pulseDutyCycle_;    // [-]
    scalar              reflectivity_;      // [-]
    scalar              transmission_;      // [-] (>=0 overrides reflectivity)
    scalar              incidenceAngle_;    // [rad]
    bool                gaussianProfile_;
    label               maxReflections_;
    bool                continuousLaser_;

    scalar              laserStartTime_;    // [s]
    scalar              laserEndTime_;      // [s]

    // Donor film bounds (y)
    scalar              filmYMin_;
    scalar              filmYMax_;

    // Off-film scaling
    scalar              offFilmScale_;
    scalar              offFilmMax_;

    // Cached source
    mutable tmp<volScalarField> tSource_;
    mutable bool         sourceValid_;
    mutable scalar       cumulativeEnergy_;
    mutable label        lastTimeIndex_;

    // --- helpers ---
    void   calculateSource() const;
    bool   validateParameters() const;
    bool   checkPhysicalBounds() const;
    scalar calculateGaussianIntensity(const scalar R) const;
    bool   isInBeam(const point& p) const;
    bool   checkEnergyConservation() const;

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

    const dimensionedScalar& peakIntensity() const { return peakIntensity_; }
    const dimensionedScalar& pulseWidth()    const { return pulseWidth_; }
    const dimensionedScalar& wavelength()    const { return wavelength_; }
    const dimensionedScalar& spotSize()      const { return spotSize_; }
    const dimensionedScalar& pulseEnergy()   const { return pulseEnergy_; }
    const point&              focus()        const { return focus_; }
    const vector&             direction()    const { return direction_; }
    scalar laserStartTime() const { return laserStartTime_; }
    scalar laserEndTime()   const { return laserEndTime_; }
    bool   isContinuous()   const { return continuousLaser_; }

    bool valid() const;
    void write() const;
};

} // End namespace Foam
#endif

