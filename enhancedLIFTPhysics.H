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
    Foam::enhancedLIFTPhysics

Description
    Combined physics models for enhanced LIFT simulations:

    1. Phase Explosion & Spinodal Decomposition
       - Explosive boiling above 0.9×T_critical (~8550 K for Ti)
       - Enhanced mass flux (100× normal evaporation)
       - Explosive pressure contribution

    2. Plasma Formation & Ionization
       - Saha equation for ionization degree (T > 30,000 K)
       - Plasma pressure and electron density
       - Laser absorption shielding

    3. Droplet/Jet Breakup (Weber Number)
       - Weber number-based breakup criterion (We > 10)
       - Droplet size prediction
       - Breakup rate calculation

    All three models are optional and controlled via controlDict switches.

References:
    - Leveugle et al., Appl. Phys. A 79, 1643-1655 (2004)
    - Perez & Lewis, J. Appl. Phys. 103, 093105 (2008)
    - Pilch & Erdman, Int. J. Multiphase Flow 13, 741-757 (1987)
    - Feinaeugle et al., Appl. Surf. Sci. 418, 2017

SourceFiles
    enhancedLIFTPhysics.C

\*---------------------------------------------------------------------------*/

#ifndef enhancedLIFTPhysics_H
#define enhancedLIFTPhysics_H

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
                    Class enhancedLIFTPhysics Declaration
\*---------------------------------------------------------------------------*/

class enhancedLIFTPhysics
{
public:

    // Public data structures

        //- Phase explosion data
        struct PhaseExplosionData
        {
            dimensionedScalar T_critical_;
            dimensionedScalar T_spinodal_;
            dimensionedScalar p_critical_;
            dimensionedScalar tau_explosion_;
            scalar explosionMultiplier_;
            Switch enabled_;

            volScalarField explosiveMassSource_;
            volScalarField explosivePressure_;
            volScalarField explosionIndicator_;

            PhaseExplosionData
            (
                const fvMesh& mesh,
                const dictionary& dict
            );
        };

        //- Plasma ionization data
        struct PlasmaData
        {
            dimensionedScalar ionizationEnergy_;
            dimensionedScalar T_ionization_;
            dimensionedScalar k_B_;
            dimensionedScalar m_e_;
            dimensionedScalar m_atom_;
            dimensionedScalar n_atom_;
            Switch enabled_;

            volScalarField ionizationDegree_;
            volScalarField n_electron_;
            volScalarField plasmaPressure_;
            volScalarField plasmaShielding_;

            PlasmaData
            (
                const fvMesh& mesh,
                const dictionary& dict
            );
        };

        //- Droplet breakup data
        struct BreakupData
        {
            scalar We_critical_;
            dimensionedScalar minJetDiameter_;
            dimensionedScalar maxJetDiameter_;
            scalar breakupTimeCoeff_;
            Switch enabled_;

            volScalarField WeberNumber_;
            volScalarField breakupIndicator_;
            volScalarField dropletDiameter_;
            volScalarField breakupRate_;

            BreakupData
            (
                const fvMesh& mesh,
                const dictionary& dict
            );
        };


private:

    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to mixture
        const twoPhaseMixtureThermo& mixture_;

        //- Phase explosion data
        autoPtr<PhaseExplosionData> phaseExplosion_;

        //- Plasma ionization data
        autoPtr<PlasmaData> plasma_;

        //- Droplet breakup data
        autoPtr<BreakupData> breakup_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        enhancedLIFTPhysics(const enhancedLIFTPhysics&);

        //- Disallow default bitwise assignment
        void operator=(const enhancedLIFTPhysics&);

        //- Calculate Saha ionization degree
        scalar sahaIonization
        (
            const scalar T,
            const scalar n_atoms
        ) const;

        //- Calculate characteristic length for breakup
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
            const volScalarField& L_char,
            const volScalarField& alpha1
        );


public:

    //- Runtime type information
    TypeName("enhancedLIFTPhysics");


    // Constructors

        //- Construct from mesh, mixture and dictionary
        enhancedLIFTPhysics
        (
            const fvMesh& mesh,
            const twoPhaseMixtureThermo& mixture,
            const dictionary& dict
        );


    //- Destructor
    virtual ~enhancedLIFTPhysics();


    // Member Functions

        // Access - Phase Explosion

            //- Is phase explosion enabled?
            bool phaseExplosionEnabled() const
            {
                return phaseExplosion_.valid() && phaseExplosion_->enabled_;
            }

            //- Return explosive mass source
            const volScalarField& explosiveMassSource() const
            {
                return phaseExplosion_->explosiveMassSource_;
            }

            //- Return explosive pressure
            const volScalarField& explosivePressure() const
            {
                return phaseExplosion_->explosivePressure_;
            }

            //- Return explosion indicator
            const volScalarField& explosionIndicator() const
            {
                return phaseExplosion_->explosionIndicator_;
            }

            //- Calculate enhanced mass flux multiplier
            scalar explosionFluxMultiplier(const scalar T) const;


        // Access - Plasma Ionization

            //- Is plasma model enabled?
            bool plasmaEnabled() const
            {
                return plasma_.valid() && plasma_->enabled_;
            }

            //- Return ionization degree
            const volScalarField& ionizationDegree() const
            {
                return plasma_->ionizationDegree_;
            }

            //- Return electron density
            const volScalarField& electronDensity() const
            {
                return plasma_->n_electron_;
            }

            //- Return plasma pressure
            const volScalarField& plasmaPressure() const
            {
                return plasma_->plasmaPressure_;
            }

            //- Return plasma shielding coefficient
            const volScalarField& plasmaShielding() const
            {
                return plasma_->plasmaShielding_;
            }


        // Access - Droplet Breakup

            //- Is breakup model enabled?
            bool breakupEnabled() const
            {
                return breakup_.valid() && breakup_->enabled_;
            }

            //- Return Weber number
            const volScalarField& WeberNumber() const
            {
                return breakup_->WeberNumber_;
            }

            //- Return breakup indicator
            const volScalarField& breakupIndicator() const
            {
                return breakup_->breakupIndicator_;
            }

            //- Return droplet diameter
            const volScalarField& dropletDiameter() const
            {
                return breakup_->dropletDiameter_;
            }


        // Evaluation

            //- Update all enabled physics models
            void updateAll
            (
                const volScalarField& T,
                const volScalarField& alpha1,
                const volScalarField& rho,
                const volVectorField& U
            );

            //- Update phase explosion model
            void updatePhaseExplosion
            (
                const volScalarField& T,
                const volScalarField& alpha1,
                const volScalarField& rho
            );

            //- Update plasma ionization model
            void updatePlasma
            (
                const volScalarField& T,
                const volScalarField& rho,
                const volScalarField& alpha1
            );

            //- Update droplet breakup model
            void updateBreakup
            (
                const volScalarField& alpha1,
                const volVectorField& U,
                const volScalarField& rho,
                const volScalarField& T
            );

            //- Apply droplet breakup to alpha field
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
