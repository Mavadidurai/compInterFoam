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
#include "fvc.H"

Foam::compressibleInterPhaseTransportModel::compressibleInterPhaseTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& rhoPhi,
    const surfaceScalarField& alphaPhi10,
    twoPhaseMixtureThermo& mixture
)
:
    twoPhaseTransport_(false),
    mixture_(mixture),
    phi_(phi),
    rhoPhi_(rhoPhi),
    alphaPhi10_(alphaPhi10)
{
    word simulationType("laminar");

    const auto readSimulationType =
        [&](const word& dictName) -> bool
        {
            IOobject io
            (
                dictName,
                U.time().constant(),
                U.db(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );

            if (!io.typeHeaderOk<IOdictionary>(true))
            {
                return false;
            }

            IOdictionary dict(io);

            if (dict.found("simulationType"))
            {
                dict.lookup("simulationType") >> simulationType;
                return true;
            }

            return false;
        };

    if (!readSimulationType(compressible::turbulenceModel::typeName))
    {
        readSimulationType("turbulenceProperties");
    }

    twoPhaseTransport_ = (simulationType == "twoPhaseTransport");

    if (twoPhaseTransport_)
    {
        const volScalarField& alpha1 = mixture_.alpha1();
        const volScalarField& alpha2 = mixture_.alpha2();

        const volScalarField& rho1 = mixture_.thermo1().rho();
        const volScalarField& rho2 = mixture_.thermo2().rho();

        alphaRhoPhi1_ = tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRhoPhi", alpha1.group()),
                    U.time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(rho1)*alphaPhi10_
            )
        );

        alphaRhoPhi2_ = tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRhoPhi", alpha2.group()),
                    U.time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(rho2)*(phi_ - alphaPhi10_)
            )
        );

turbulence1_ = PhaseCompressibleTurbulenceModel<rhoThermo>::New
        (
            alpha1,
            rho1,
            U,
            alphaRhoPhi1_(),
            phi_,
            mixture.thermo1()
        );

turbulence2_ = PhaseCompressibleTurbulenceModel<rhoThermo>::New
        (
            alpha2,
            rho2,
            U,
            alphaRhoPhi2_(),
            phi_,
            mixture.thermo2()
        );
    }
    else
    {
        turbulence_ = compressible::turbulenceModel::New
        (
            rho,
            U,
            rhoPhi_,
            mixture
        );
    }
}

Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseTransportModel::alphaEff() const
{
    if (twoPhaseTransport_)
    {
        return
            mixture_.alpha1()*mixture_.thermo1().alphaEff
            (
                turbulence1_->alphatEff()
            )
          + mixture_.alpha2()*mixture_.thermo2().alphaEff
            (
                turbulence2_->alphatEff()
            );
    }

    return mixture_.alphaEff(turbulence_->alphat());
}

Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseTransportModel::kappaEff() const
{
    if (twoPhaseTransport_)
    {
        return
            mixture_.alpha1()*mixture_.thermo1().kappaEff
            (
                turbulence1_->alphatEff()
            )
          + mixture_.alpha2()*mixture_.thermo2().kappaEff
            (
                turbulence2_->alphatEff()
            );
    }

    return mixture_.kappaEff(turbulence_->alphat());
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::compressibleInterPhaseTransportModel::divDevRhoReff
(
    volVectorField& U
) const
{
    if (twoPhaseTransport_)
    {
        return turbulence1_->divDevRhoReff(U) + turbulence2_->divDevRhoReff(U);
    }

    return turbulence_->divDevRhoReff(U);
}

void Foam::compressibleInterPhaseTransportModel::correctPhasePhi()
{
    if (!twoPhaseTransport_)
    {
        return;
    }

    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    alphaRhoPhi1_.ref() = fvc::interpolate(rho1)*alphaPhi10_;
    alphaRhoPhi2_.ref() = fvc::interpolate(rho2)*(phi_ - alphaPhi10_);
}

void Foam::compressibleInterPhaseTransportModel::correct()
{
    if (twoPhaseTransport_)
    {
        turbulence1_->correct();
        turbulence2_->correct();
    }
    else
    {
        turbulence_->correct();
    }
}
// ************************************************************************* //
