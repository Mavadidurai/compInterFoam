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
    Foam::plasmaIonizationModel

Description
    Ionization and plasma pressure model for ultrafast laser-induced plasmas.

    Above the ionization threshold (~few eV, ~30,000-50,000 K), electrons
    are liberated from atoms creating a plasma. This model computes:
    - Ionization degree (Saha equation)
    - Plasma pressure (electron + ion contributions)
    - Plasma shielding effects on laser absorption

    The plasma formation is important for:
    - Pressures above the normal recoil pressure
    - Shielding/absorption of subsequent laser pulses
    - Enhanced material removal rates

    References:
    - Perez & Lewis, J. Appl. Phys. 103, 093105 (2008)
    - Bulgakova et al., Phys. Rev. E 62, 5624 (2000)

SourceFiles
    plasmaIonizationModel.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaIonizationModel_H
#define plasmaIonizationModel_H

#include "volFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "fvMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaIonizationModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaIonizationModel
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Ionization energy [J] (eV converted to J)
        dimensionedScalar ionizationEnergy_;

        //- Ionization threshold temperature [K]
        dimensionedScalar T_ionization_;

        //- Boltzmann constant [J/K]
        dimensionedScalar k_B_;

        //- Electron mass [kg]
        dimensionedScalar m_e_;

        //- Atomic mass [kg]
        dimensionedScalar m_atom_;

        //- Number density of atoms [1/m^3]
        dimensionedScalar n_atom_;

        //- Enable/disable plasma model
        Switch enabled_;

        //- Ionization degree field (0-1)
        volScalarField ionizationDegree_;

        //- Electron number density [1/m^3]
        volScalarField n_electron_;

        //- Plasma pressure [Pa]
        volScalarField plasmaPressure_;

        //- Plasma shielding coefficient (0-1, reduces laser absorption)
        volScalarField plasmaShielding_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        plasmaIonizationModel(const plasmaIonizationModel&);

        //- Disallow default bitwise assignment
        void operator=(const plasmaIonizationModel&);

        //- Calculate ionization degree using Saha equation
        scalar sahaIonization
        (
            const scalar T,
            const scalar n_atoms
        ) const;


public:

    //- Runtime type information
    TypeName("plasmaIonizationModel");


    // Constructors

        //- Construct from mesh and dictionary
        plasmaIonizationModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~plasmaIonizationModel();


    // Member Functions

        // Access

            //- Return ionization degree
            const volScalarField& ionizationDegree() const
            {
                return ionizationDegree_;
            }

            //- Return electron density
            const volScalarField& electronDensity() const
            {
                return n_electron_;
            }

            //- Return plasma pressure
            const volScalarField& plasmaPressure() const
            {
                return plasmaPressure_;
            }

            //- Return plasma shielding coefficient
            const volScalarField& plasmaShielding() const
            {
                return plasmaShielding_;
            }

            //- Is plasma model enabled?
            bool enabled() const
            {
                return enabled_;
            }


        // Evaluation

            //- Update plasma fields based on current temperature and density
            void update
            (
                const volScalarField& T,
                const volScalarField& rho,
                const volScalarField& alpha1
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
