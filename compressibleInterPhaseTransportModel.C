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
        // Use nut() instead of alphat() for OpenFOAM v2406 compatibility
        const tmp<volScalarField> nut1Tmp = turbulence1_->nut();
        const volScalarField& nut1 = nut1Tmp();

        const tmp<volScalarField> nut2Tmp = turbulence2_->nut();
        const volScalarField& nut2 = nut2Tmp();

        // Convert kinematic turbulent viscosity to thermal diffusivity
        // using Prandtl number approach: alphat = nut/Prt
        const scalar Prt1 = 0.85; // Typical turbulent Prandtl number for phase 1
        const scalar Prt2 = 0.85; // Typical turbulent Prandtl number for phase 2
        
        const tmp<volScalarField> alphat1 = nut1/Prt1;
        const tmp<volScalarField> alphat2 = nut2/Prt2;

        return
            mixture_.alpha1()*mixture_.thermo1().alphaEff(alphat1())
          + mixture_.alpha2()*mixture_.thermo2().alphaEff(alphat2());
    }

    const tmp<volScalarField> nutTmp = turbulence_->nut();
    const volScalarField& nut = nutTmp();
    
    // Convert to thermal diffusivity using Prandtl number
    const scalar Prt = 0.85; // Typical turbulent Prandtl number
    const tmp<volScalarField> alphat = nut/Prt;
    
    return mixture_.alphaEff(alphat());
}

Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseTransportModel::kappaEff() const
{
    if (twoPhaseTransport_)
    {
        // Use nut() instead of alphat() for OpenFOAM v2406 compatibility
        const tmp<volScalarField> nut1Tmp = turbulence1_->nut();
        const volScalarField& nut1 = nut1Tmp();

        const tmp<volScalarField> nut2Tmp = turbulence2_->nut();
        const volScalarField& nut2 = nut2Tmp();

        // Convert kinematic turbulent viscosity to thermal diffusivity
        const scalar Prt1 = 0.85; // Typical turbulent Prandtl number for phase 1
        const scalar Prt2 = 0.85; // Typical turbulent Prandtl number for phase 2
        
        const tmp<volScalarField> alphat1 = nut1/Prt1;
        const tmp<volScalarField> alphat2 = nut2/Prt2;

        return
            mixture_.alpha1()*mixture_.thermo1().kappaEff(alphat1())
          + mixture_.alpha2()*mixture_.thermo2().kappaEff(alphat2());
    }

    const tmp<volScalarField> nutTmp = turbulence_->nut();
    const volScalarField& nut = nutTmp();
    
    // Convert to thermal diffusivity using Prandtl number
    const scalar Prt = 0.85; // Typical turbulent Prandtl number
    const tmp<volScalarField> alphat = nut/Prt;
    
    return mixture_.kappaEff(alphat());
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
