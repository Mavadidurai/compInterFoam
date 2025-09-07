/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::compressibleInterPhaseTransportModel
Description
    Transport model selection class for the compressibleInterFoam family of
    solvers.
    By default the standard mixture transport modelling approach is used in
    which a single momentum stress model (laminar, non-Newtonian, LES or RAS) is
    constructed for the mixture.  However if the \c simulationType in
    constant/turbulenceProperties is set to \c twoPhaseTransport the alternative
    Euler-Euler two-phase transport modelling approach is used in which separate
    stress models (laminar, non-Newtonian, LES or RAS) are instantiated for each
    of the two phases allowing for different modeling for the phases.
SourceFiles
    compressibleInterPhaseTransportModel.C
\*---------------------------------------------------------------------------*/
#ifndef compressibleInterPhaseTransportModel_H
#define compressibleInterPhaseTransportModel_H
#include "twoPhaseMixture.H"
#include "twoPhaseMixtureThermo.H"
#include "turbulentFluidThermoModel.H"    
#include "compressibleTurbulenceModel.H"   
#include "fluidThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
             Class compressibleInterPhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/
class compressibleInterPhaseTransportModel
{
    // Private data
    const twoPhaseMixtureThermo& mixture_;
    const surfaceScalarField& phi_;
    const surfaceScalarField& rhoPhi_;
    const surfaceScalarField& alphaPhi10_;
    // Single turbulence model for the mixture
    autoPtr<compressible::turbulenceModel> turbulence_;
public:
    // Constructors
    compressibleInterPhaseTransportModel
    (
        const volScalarField& rho,
        const volVectorField& U,
        const surfaceScalarField& phi,
        const surfaceScalarField& rhoPhi,
        const surfaceScalarField& alphaPhi10,
        const twoPhaseMixtureThermo& mixture
    );
    // Member Functions
    tmp<volScalarField> alphaEff() const;
    tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;
    void correctPhasePhi() {}  // Empty for simplified model
    void correct();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
