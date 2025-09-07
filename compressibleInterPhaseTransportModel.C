/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
\*---------------------------------------------------------------------------*/
#include "compressibleInterPhaseTransportModel.H"
Foam::compressibleInterPhaseTransportModel::compressibleInterPhaseTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& rhoPhi,
    const surfaceScalarField& alphaPhi10,
    const twoPhaseMixtureThermo& mixture
)
:
    mixture_(mixture),
    phi_(phi),
    rhoPhi_(rhoPhi), 
    alphaPhi10_(alphaPhi10),
    turbulence_
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            rhoPhi,
            mixture
        )
    )
{
    turbulence_->validate();
}
Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseTransportModel::alphaEff() const
{
    return mixture_.alphaEff(turbulence_->alphat());
}
Foam::tmp<Foam::fvVectorMatrix>
Foam::compressibleInterPhaseTransportModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return turbulence_->divDevRhoReff(U);
}
void Foam::compressibleInterPhaseTransportModel::correct()
{
    turbulence_->correct();
}
// ************************************************************************* //
