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
    Foam::phaseExplosionModel

Description
    Phase explosion and spinodal decomposition model for ultrafast laser
    heating beyond the Clausius-Clapeyron regime.

    When the temperature exceeds ~0.9*T_critical, the metastable liquid
    undergoes explosive boiling (phase explosion). This model computes:
    - Spinodal limit temperature
    - Enhanced mass flux above spinodal threshold
    - Explosive pressure rise

    References:
    - Leveugle et al., Appl. Phys. A 79, 1643-1655 (2004)
    - Zhigilei et al., Appl. Phys. A 114, 11-15 (2014)

SourceFiles
    phaseExplosionModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseExplosionModel_H
#define phaseExplosionModel_H

#include "volFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class phaseExplosionModel Declaration
\*---------------------------------------------------------------------------*/

class phaseExplosionModel
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Critical temperature [K]
        dimensionedScalar T_critical_;

        //- Spinodal temperature (0.9 * T_critical) [K]
        dimensionedScalar T_spinodal_;

        //- Critical pressure [Pa]
        dimensionedScalar p_critical_;

        //- Phase explosion time scale [s]
        dimensionedScalar tau_explosion_;

        //- Enhanced mass flux multiplier (dimensionless)
        scalar explosionMultiplier_;

        //- Enable/disable phase explosion model
        Switch enabled_;

        //- Explosive mass source field [kg/m^3/s]
        volScalarField explosiveMassSource_;

        //- Explosive pressure field [Pa]
        volScalarField explosivePressure_;

        //- Phase explosion indicator (0 = normal, 1 = explosive)
        volScalarField explosionIndicator_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseExplosionModel(const phaseExplosionModel&);

        //- Disallow default bitwise assignment
        void operator=(const phaseExplosionModel&);


public:

    //- Runtime type information
    TypeName("phaseExplosionModel");


    // Constructors

        //- Construct from mesh and dictionary
        phaseExplosionModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseExplosionModel();


    // Member Functions

        // Access

            //- Return critical temperature
            const dimensionedScalar& T_critical() const
            {
                return T_critical_;
            }

            //- Return spinodal temperature
            const dimensionedScalar& T_spinodal() const
            {
                return T_spinodal_;
            }

            //- Return explosive mass source
            const volScalarField& explosiveMassSource() const
            {
                return explosiveMassSource_;
            }

            //- Return explosive pressure
            const volScalarField& explosivePressure() const
            {
                return explosivePressure_;
            }

            //- Return explosion indicator
            const volScalarField& explosionIndicator() const
            {
                return explosionIndicator_;
            }

            //- Is phase explosion enabled?
            bool enabled() const
            {
                return enabled_;
            }


        // Evaluation

            //- Update phase explosion fields based on current temperature
            void update
            (
                const volScalarField& T,
                const volScalarField& alpha1,
                const volScalarField& rho
            );

            //- Calculate enhanced mass flux multiplier
            scalar explosionFluxMultiplier(const scalar T) const;

            //- Calculate explosive pressure
            scalar explosionPressure
            (
                const scalar T,
                const scalar rho
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
