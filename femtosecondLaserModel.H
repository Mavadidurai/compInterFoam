/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom LIFT Solver
    
    Description
    Header file for femtosecond laser model used in LIFT simulation.
    Provides laser-material interaction with substrate-backed geometry.

\*---------------------------------------------------------------------------*/

#ifndef femtosecondLaserModel_H
#define femtosecondLaserModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "mathematicalConstants.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class femtosecondLaserModel Declaration
\*---------------------------------------------------------------------------*/

class femtosecondLaserModel
{
private:

    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;
        
        //- Reference to dictionary
        const dictionary& dict_;

        // Laser parameters
        dimensionedScalar peakIntensity_;    // Peak intensity [W/m²]
        dimensionedScalar pulseWidth_;       // Pulse duration [s]
        dimensionedScalar wavelength_;       // Laser wavelength [m]
        dimensionedScalar absorptionCoeff_;  // Absorption coefficient [1/m]
        dimensionedScalar spotSize_;         // Laser spot size [m]
        dimensionedScalar pulseEnergy_;      // Pulse energy [J]
        
        // Geometrical parameters
        vector direction_;              // Laser propagation direction
        point focus_;                   // Focal point coordinates
        point initialFocus_;            // Stored initial focus for scanning
        vector scanVelocity_;           // Beam scanning velocity [m/s]
        scalar pulseFrequency_;         // Pulse repetition frequency [Hz]
        scalar pulseDutyCycle_;         // Fraction of period laser is on
        
        // Optional parameters
        scalar reflectivity_;           // Surface reflectivity [-]
        bool gaussianProfile_;          // Use Gaussian spatial profile
        label maxReflections_;          // Maximum number of reflections
        
        // Temporal parameters
        bool continuousLaser_;          // Continuous vs pulsed mode
        scalar laserStartTime_;         // Laser start time [s]
        scalar laserEndTime_;           // Laser end time [s]
        //- Donor film bounds in y-direction [m]
        scalar filmYMin_;              // Lower y bound of metal film
        scalar filmYMax_;              // Upper y bound of metal film
        //- Scaling factor and maximum intensity when outside donor film
        scalar offFilmScale_;          // Intensity scaling outside film [-]
        scalar offFilmMax_;            // Max volumetric intensity outside film [W/m^3]

        // Cached fields
        mutable tmp<volScalarField> tSource_;  // Laser source field
        mutable bool sourceValid_;             // Source validity flag
        
        //- Running total of energy deposited [J]
        mutable scalar cumulativeEnergy_;

        //- Last processed time index for restart detection
        mutable label lastTimeIndex_;

     // Private Member Functions

        //- Calculate laser source term
        void calculateSource() const;

        //- Validate laser parameters
        bool validateParameters() const;

        //- Check physical bounds of parameters
        bool checkPhysicalBounds() const;

        //- Calculate Gaussian intensity profile
        scalar calculateGaussianIntensity
        (
            const scalar R,      // Radial distance from beam center
            const scalar z       // Axial distance along beam
        ) const;

        //- Check if point is within laser beam
        bool isInBeam(const point& p) const;
        
        //- Check energy conservation
        bool checkEnergyConservation() const;

        //- Disallow default bitwise copy construct
        femtosecondLaserModel(const femtosecondLaserModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const femtosecondLaserModel&) = delete;

public:

    //- Runtime type information
    TypeName("femtosecondLaserModel");

    // Constructors

        //- Construct from mesh and dictionary
        femtosecondLaserModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~femtosecondLaserModel() = default;

    // Member Functions

        // Evolution

            //- Update the laser model (invalidate cached data)
            void update();

            //- Correct the laser model (calculate source)
            virtual void correct()
            {
                sourceValid_ = false;
                calculateSource();
            }
            //- Correct the laser model for a given current time, allowing
            //  time-dependent beam motion or intensity modulation
            virtual void correct(const scalar currentTime);

        // Access

            //- Return the laser source term [W/m³]
            tmp<volScalarField> source() const;

            //- Return peak intensity
            const dimensionedScalar& peakIntensity() const
            {
                return peakIntensity_;
            }

            //- Return pulse width
            const dimensionedScalar& pulseWidth() const
            {
                return pulseWidth_;
            }

            //- Return wavelength
            const dimensionedScalar& wavelength() const
            {
                return wavelength_;
            }

            //- Return spot size
            const dimensionedScalar& spotSize() const
            {
                return spotSize_;
            }

            //- Return pulse energy
            const dimensionedScalar& pulseEnergy() const
            {
                return pulseEnergy_;
            }

            //- Return focal point
            const point& focus() const
            {
                return focus_;
            }

            //- Return laser direction
            const vector& direction() const
            {
                return direction_;
            }

            //- Return start time
            scalar laserStartTime() const
            {
                return laserStartTime_;
            }

            //- Return end time
            scalar laserEndTime() const
            {
                return laserEndTime_;
            }

            //- Check if laser is continuous
            bool isContinuous() const
            {
                return continuousLaser_;
            }

        // Validation

            //- Check if model parameters are valid
            bool valid() const;

        // I-O

            //- Write laser model information
            void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
