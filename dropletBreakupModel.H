/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Class
    Foam::dropletBreakupModel

Description
    Weber-number based droplet/jet breakup model for LIFT simulations.

    When molten jets or ligaments form, they may fragment into smaller
    droplets based on the Weber number:
        We = ρ * U^2 * L / σ

    where:
        ρ = density
        U = relative velocity
        L = characteristic length (jet/droplet diameter)
        σ = surface tension

    Breakup occurs when We > We_critical (typically ~10-12).

    This model:
    - Detects high-velocity jets/ligaments
    - Calculates local Weber number
    - Applies interface smoothing/breakup where We > We_crit
    - Tracks breakup events and droplet sizes

    References:
    - Pilch & Erdman, Int. J. Multiphase Flow 13, 741-757 (1987)
    - Guildenbecher et al., Int. J. Multiphase Flow 35, 741-757 (2009)

SourceFiles
    dropletBreakupModel.C

\*---------------------------------------------------------------------------*/

#ifndef dropletBreakupModel_H
#define dropletBreakupModel_H

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "fvMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class twoPhaseMixtureThermo;

/*---------------------------------------------------------------------------*\
                    Class dropletBreakupModel Declaration
\*---------------------------------------------------------------------------*/

class dropletBreakupModel
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to mixture
        const twoPhaseMixtureThermo& mixture_;

        //- Critical Weber number for breakup
        scalar We_critical_;

        //- Minimum jet diameter for breakup consideration [m]
        dimensionedScalar minJetDiameter_;

        //- Maximum jet diameter for breakup consideration [m]
        dimensionedScalar maxJetDiameter_;

        //- Breakup timescale coefficient
        scalar breakupTimeCoeff_;

        //- Enable/disable breakup model
        Switch enabled_;

        //- Local Weber number field
        volScalarField WeberNumber_;

        //- Breakup indicator (0 = stable, 1 = breaking)
        volScalarField breakupIndicator_;

        //- Estimated droplet diameter after breakup [m]
        volScalarField dropletDiameter_;

        //- Breakup rate [1/s]
        volScalarField breakupRate_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dropletBreakupModel(const dropletBreakupModel&);

        //- Disallow default bitwise assignment
        void operator=(const dropletBreakupModel&);

        //- Calculate characteristic length scale (jet diameter)
        tmp<volScalarField> calculateCharacteristicLength
        (
            const volScalarField& alpha1
        ) const;

        //- Calculate Weber number
        void calculateWeberNumber
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& sigma,
            const volScalarField& L_char
        );


public:

    //- Runtime type information
    TypeName("dropletBreakupModel");


    // Constructors

        //- Construct from components
        dropletBreakupModel
        (
            const fvMesh& mesh,
            const twoPhaseMixtureThermo& mixture,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dropletBreakupModel();


    // Member Functions

        // Access

            //- Return Weber number field
            const volScalarField& WeberNumber() const
            {
                return WeberNumber_;
            }

            //- Return breakup indicator
            const volScalarField& breakupIndicator() const
            {
                return breakupIndicator_;
            }

            //- Return droplet diameter
            const volScalarField& dropletDiameter() const
            {
                return dropletDiameter_;
            }

            //- Return breakup rate
            const volScalarField& breakupRate() const
            {
                return breakupRate_;
            }

            //- Is breakup model enabled?
            bool enabled() const
            {
                return enabled_;
            }


        // Evaluation

            //- Update breakup fields
            void update
            (
                const volScalarField& alpha1,
                const volVectorField& U,
                const volScalarField& rho
            );

            //- Apply breakup to alpha field (interface remapping)
            void applyBreakup
            (
                volScalarField& alpha1,
                const volVectorField& U
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
