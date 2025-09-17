/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "twoPhaseMixtureThermo.H"
#include "interfaceProperties.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "collatedFileOperation.H"
#include "dimensionSets.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseMixtureThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureThermo::twoPhaseMixtureThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    psiThermo(U.mesh(), word::null),
    twoPhaseMixture(U.mesh(), *this),
    interfaceProperties(alpha1(), U, *this),
    thermo1_(nullptr),
    thermo2_(nullptr),
   // Initialise phase-change properties; values populated from metal dictionary
    latentHeat_(0.0),
    T_melt_(0.0),
    T_vapor_(0.0),
    dtFloor_(1e-12),
    Q_laser_
    (
        IOobject
        (
            "Q_laser",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("Q0", dimPower/dimVolume, 0.0)  // [W/m³]
    ),
    phaseChangeSource_
    (
        IOobject
        (
            "phaseChangeSource",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("source", dimTemperature/dimTime, 0.0)  // [K/s]
       ),
     phaseChangeRelaxCoeff_
    (
        IOobject
        (
            "phaseChangeRelaxCoeff",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("relax", dimless/dimTime, 0.0)
    ),
    dgdt_
    (
        IOobject
        (
            "dgdt",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("dgdt", dimless/dimTime, 0.0)
    ),
    nu1_("nu1", dimViscosity, 0.0),
    nu2_("nu2", dimViscosity, 0.0),
    rho1_("rho1", dimDensity, 0.0),
    rho2_("rho2", dimDensity, 0.0)
{
        // Read metal phase thermophysical properties from dedicated dictionary
    IOdictionary metalDict
    (
        IOobject
        (
            "thermophysicalProperties.metal",
            U.mesh().time().constant(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    latentHeat_ = metalDict.lookupOrDefault<scalar>("hf", 435e3);
    T_melt_ = metalDict.lookupOrDefault<scalar>("Tsol", 1941.0);
    T_vapor_ = metalDict.lookupOrDefault<scalar>("Tvap", 3000.0);

    {
        volScalarField T1(IOobject::groupName("T", phase1Name()), T_);
        T1.write();
    }

    {
        volScalarField T2(IOobject::groupName("T", phase2Name()), T_);
        T2.write();
    }

    // Note: we're writing files to be read in immediately afterwards.
    //       Avoid any thread-writing problems.
    fileHandler().flush();
        const dictionary& transportDict =
        U.mesh().lookupObject<dictionary>("transportProperties");

    const dictionary& phase1Dict = transportDict.subDict(phase1Name());
    const dictionary& phase2Dict = transportDict.subDict(phase2Name());

    nu1_ = phase1Dict.lookupOrDefault<dimensionedScalar>
    (
        "nu",
        dimensionedScalar("nu", dimViscosity, 0)
    );
    nu2_ = phase2Dict.lookupOrDefault<dimensionedScalar>
    (
        "nu",
        dimensionedScalar("nu", dimViscosity, 0)
    );

    rho1_ = phase1Dict.lookupOrDefault<dimensionedScalar>
    (
        "rho",
        dimensionedScalar("rho", dimDensity, 0)
    );
    rho2_ = phase2Dict.lookupOrDefault<dimensionedScalar>
    (
        "rho",
        dimensionedScalar("rho", dimDensity, 0)
    );

    Info<< "Transport properties:" << nl
        << "    " << phase1Name() << ": nu=" << nu1_.value()
        << ", rho=" << rho1_.value() << nl
        << "    " << phase2Name() << ": nu=" << nu2_.value()
        << ", rho=" << rho2_.value() << endl;


    thermo1_ = rhoThermo::New(U.mesh(), phase1Name());
    thermo2_ = rhoThermo::New(U.mesh(), phase2Name());

    correct();
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureThermo::~twoPhaseMixtureThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::twoPhaseMixtureThermo::correctThermo()
{
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();
}
void Foam::twoPhaseMixtureThermo::correct()
{
    psi_ = alpha1()*thermo1_->psi() + alpha2()*thermo2_->psi();
    mu_ = alpha1()*thermo1_->mu() + alpha2()*thermo2_->mu();
    alpha_ = alpha1()*thermo1_->alpha() + alpha2()*thermo2_->alpha();

    interfaceProperties::correct();
    // Laser heating is now supplied externally via setQLaser()

    phaseChangeSource_ = computePhaseChange()();
    phaseChangeSource_.correctBoundaryConditions();
    phaseChangeRelaxCoeff_.correctBoundaryConditions();

    // Compute mass-transfer rate [1/s]
    dgdt_ = computeMassTransfer()();
    dgdt_.correctBoundaryConditions();

}
void Foam::twoPhaseMixtureThermo::setQLaser(const volScalarField& src)
{
    Q_laser_ = src;
    Q_laser_.correctBoundaryConditions();
}
Foam::word Foam::twoPhaseMixtureThermo::thermoName() const
{
    return thermo1_->thermoName() + ',' + thermo2_->thermoName();
}
bool Foam::twoPhaseMixtureThermo::incompressible() const
{
    return thermo1_->incompressible() && thermo2_->incompressible();
}
Foam::dimensionedScalar Foam::twoPhaseMixtureThermo::latentHeat() const
{
    return Foam::dimensionedScalar
    (
        "latentHeat",
        Foam::dimLatentHeat,
        latentHeat_
    );
}
Foam::tmp<Foam::volScalarField>
Foam::twoPhaseMixtureThermo::sigma() const
{
    return this->Foam::interfaceProperties::sigmaK();
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::computePhaseChange()
{
    tmp<volScalarField> tSource
    (
        new volScalarField
        (
            IOobject
            (
                "phaseChangeSource",
                T_.time().timeName(),
                T_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionedScalar("source", dimTemperature/dimTime, 0.0)
        )
    );

    volScalarField& source = tSource.ref();
    // Reset the implicit relaxation coefficient
    phaseChangeRelaxCoeff_ = dimensionedScalar("relax", dimless/dimTime, 0.0);
    // Access coefficients from transportProperties
    const dictionary& transportDict =
        T_.mesh().lookupObject<dictionary>("transportProperties");

    if (!transportDict.found("phaseChangeCoeffs"))
    {
        FatalErrorInFunction
            << "phaseChangeCoeffs not found in transportProperties" << nl
            << "Supply a phaseChangeCoeffs sub-dictionary or disable phase change." << nl
            << exit(FatalError);
        return tSource;
    }

    const dictionary& pc = transportDict.subDict("phaseChangeCoeffs");

    const scalar Tvapor = pc.lookupOrDefault<scalar>("Tvapor", T_vapor_);
    const scalar windowWidth = pc.lookupOrDefault<scalar>("windowWidth", 0.0);
    dtFloor_ = pc.lookupOrDefault<scalar>("dtFloor", dtFloor_);
    scalar relaxationRate = pc.lookupOrDefault<scalar>("relaxationRate", -1.0);
    if (relaxationRate < 0)
    {
        const scalar relaxationTime =
            pc.lookupOrDefault<scalar>("relaxationTime", -1.0);

        if (relaxationTime > 0)
        {
            relaxationRate = 1.0/relaxationTime;
        }
    }

    if (relaxationRate < 0)
    {
        FatalErrorInFunction
            << "phaseChangeCoeffs requires either 'relaxationRate' [1/s]"
            << " or 'relaxationTime' [s] to be specified" << nl
            << exit(FatalError);
    }

    scalar maxSource = pc.lookupOrDefault<scalar>("maxSource", GREAT);
    if (maxSource >= GREAT && pc.found("minCoefficient"))
    {
        maxSource = pc.lookup<scalar>("minCoefficient");
    }
    const Switch onlyAboveVapor
    (
        pc.lookupOrDefault<Switch>("onlyAboveVapor", false)
    );

    List<Tuple2<scalar, scalar>> actTimes;
    if (pc.found("activationTime"))
    {
        pc.lookup("activationTime") >> actTimes;
    }
    // Log the coefficients only once
    static bool loggedPhaseChange = false;
    if (!loggedPhaseChange)
    {
        Info<< "phaseChangeCoeffs:" << nl
            << "    Tvapor        " << Tvapor << nl
            << "    windowWidth   " << windowWidth << nl
            << "    dtFloor       " << dtFloor_ << nl
            << "    relaxationRate " << relaxationRate << nl
            << "    relaxationTime "
            << (relaxationRate > SMALL ? 1.0/relaxationRate : GREAT) << nl
            << "    maxSource     " << maxSource << nl
            << "    onlyAboveVapor " << onlyAboveVapor << nl;
        if (actTimes.size())
        {
            Info<< "    activationTime";
            forAll(actTimes, i)
            {
                Info<< " (" << actTimes[i].first() << ' ' << actTimes[i].second() << ')';
            }
            Info<< nl;
        }
        else
        {
            Info<< "    activationTime none" << nl;
        }
        loggedPhaseChange = true;
    }
    bool active = true;
    if (actTimes.size())
    {
        active = false;
        const scalar timeVal = T_.time().value();
        forAll(actTimes, i)
        {
            if
            (
                timeVal >= actTimes[i].first()
             && timeVal <= actTimes[i].second()
            )
            {
                active = true;
                break;
            }
        }
        if (!active)
        {
            Info<< "phaseChangeCoeffs inactive at time " << timeVal << nl;
            return tSource;
        }
    }

    const scalar minAlpha = 0.0;
    const scalar maxAlpha = 1.0;

    forAll(T_, cellI)
    {
        const scalar a1 = Foam::min(Foam::max(alpha1()[cellI], minAlpha), maxAlpha);
        const scalar Tcell = T_[cellI];

        scalar available = 0.0;
        if (Tcell > Tvapor)
        {
            available = a1;
        }
        else if (Tcell < Tvapor)
        {
            available = 1.0 - a1;
        }

        scalar localRate = relaxationRate*available;

        if (windowWidth > SMALL)
        {
            localRate *= Foam::min(Foam::mag(Tcell - Tvapor)/windowWidth, 1.0);
        }

        scalar sourceVal = localRate*(Tvapor - Tcell);

        if (maxSource < GREAT)
        {
            sourceVal = Foam::max(Foam::min(sourceVal, maxSource), -maxSource);

            if (Foam::mag(Tvapor - Tcell) > SMALL)
            {
                localRate = Foam::mag(sourceVal)/(Foam::mag(Tvapor - Tcell));
            }
            else
            {
                localRate = 0.0;
            }
        }

        if (onlyAboveVapor && Tcell < Tvapor)
        {
            localRate = 0.0;
            sourceVal = 0.0;
        }

        phaseChangeRelaxCoeff_[cellI] = localRate;
        source[cellI] = sourceVal;
    }

    return tSource;
}

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseMixtureThermo::computeMassTransfer() const
{
    tmp<volScalarField> tDgdt
    (
        new volScalarField
        (
            IOobject
            (
                "dgdt",
                T_.time().timeName(),
                T_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionedScalar("dgdt", dimless/dimTime, 0.0)
        )
    );

    volScalarField& dgdt = tDgdt.ref();

    const dictionary& transportDict =
        T_.mesh().lookupObject<dictionary>("transportProperties");

    if (!transportDict.found("massTransferCoeffs"))
    {
        Info<< "massTransferCoeffs not found - skipping mass transfer" << nl;
        return tDgdt;
    }

    const dictionary& mt = transportDict.subDict("massTransferCoeffs");

    const scalar rateMax = mt.lookupOrDefault<scalar>("rateMax", -1.0);

    // lattice heat capacity [J/m^3/K] from two-temperature properties
    const dictionary& twoTempDict =
        T_.mesh().time().controlDict().subDict("twoTemperatureProperties");

    const dimensionedScalar Cl
    (
        twoTempDict.lookupOrDefault<dimensionedScalar>
        (
            "Cl",
            dimensionedScalar("Cl", dimEnergy/dimVolume/dimTemperature, 0.0)
        )
    );

    const dimensionedScalar L = latentHeat();
    const scalar ClVal = Cl.value();
    const scalar LVal = L.value();

    List<scalar> tStart, tEnd;
    if (mt.found("tStart"))
    {
        mt.lookup("tStart") >> tStart;
    }
    if (mt.found("tEnd"))
    {
        mt.lookup("tEnd") >> tEnd;
    }
    static bool loggedMassTransfer = false;
    if (!loggedMassTransfer)
    {
        Info<< "massTransferCoeffs:" << nl;
        if (rateMax > 0)
        {
            Info<< "    rateMax   " << rateMax << nl;
        }
        else
        {
            Info<< "    rateMax   none" << nl;
        }
        if (tStart.size() && tEnd.size())
        {
            Info<< "    activationTime";
            const label nWin = min(tStart.size(), tEnd.size());
            for (label i = 0; i < nWin; ++i)
            {
                Info<< " (" << tStart[i] << ' ' << tEnd[i] << ')';
            }
            Info<< nl;
        }
        else
        {
            Info<< "    activationTime none" << nl;
        }
        loggedMassTransfer = true;
    }
    bool active = true;
    if (tStart.size() && tEnd.size())
    {
        active = false;
        const scalar timeVal = T_.time().value();
        const label n = min(tStart.size(), tEnd.size());
        for (label i = 0; i < n; ++i)
        {
            if (timeVal >= tStart[i] && timeVal <= tEnd[i])
            {
                active = true;
                break;
            }
        }
        if (!active)
        {
            Info<< "massTransferCoeffs inactive at time " << timeVal << nl;
            return tDgdt;
        }
    }

    // Retrieve phase-1 density field and store temporary to avoid referencing
    // destroyed objects
    const tmp<volScalarField> rho1Tmp = thermo1_->rho();
    const volScalarField& rho1Field = rho1Tmp();

    forAll(dgdt, cellI)
    {
        scalar localRate =
            -(ClVal*phaseChangeSource_[cellI])
            /(LVal*max(rho1Field[cellI], SMALL));

        if (rateMax > 0)
        {
            localRate = max(min(localRate, rateMax), -rateMax);
        }

        dgdt[cellI] = localRate;
    }

    return tDgdt;
}
bool Foam::twoPhaseMixtureThermo::isochoric() const

{
    return thermo1_->isochoric() && thermo2_->isochoric();
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->he(p, T) + alpha2()*thermo2_->he(p, T);
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->he(p, T, cells)
      + scalarField(alpha2(), cells)*thermo2_->he(p, T, cells);
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->he(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hc() const
{
    return alpha1()*thermo1_->hc() + alpha2()*thermo2_->hc();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    WarningInFunction
        << "THE() is unsupported; returning starting temperature" << endl;
    return T0;
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    WarningInFunction
        << "THE() is unsupported; returning starting temperature" << endl;
    return T0;
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::Cp() const
{
    return alpha1()*thermo1_->Cp() + alpha2()*thermo2_->Cp();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cp(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cp(p, T, patchi);
}

Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->Cp(p, T, cells)
      + scalarField(alpha2(), cells)*thermo2_->Cp(p, T, cells);
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::Cv() const
{
    return alpha1()*thermo1_->Cv() + alpha2()*thermo2_->Cv();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cv(p, T, patchi);
}

Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::rhoEoS
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->rhoEoS(p, T, cells)
      + scalarField(alpha2(), cells)*thermo2_->rhoEoS(p, T, cells);
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::gamma() const
{
    return alpha1()*thermo1_->gamma() + alpha2()*thermo2_->gamma();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->gamma(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->gamma(p, T, patchi);
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::Cpv() const
{
    return alpha1()*thermo1_->Cpv() + alpha2()*thermo2_->Cpv();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cpv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cpv(p, T, patchi);
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::CpByCpv() const
{
    return
        alpha1()*thermo1_->CpByCpv()
      + alpha2()*thermo2_->CpByCpv();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->CpByCpv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->CpByCpv(p, T, patchi);
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::W() const
{
    return alpha1()*thermo1_->W() + alpha2()*thermo2_->W();
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::nu() const
{
    return mu()/(alpha1()*thermo1_->rho() + alpha2()*thermo2_->rho());
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::nu
(
    const label patchi
) const
{
    return
        mu(patchi)
       /(
            alpha1().boundaryField()[patchi]*thermo1_->rho(patchi)
          + alpha2().boundaryField()[patchi]*thermo2_->rho(patchi)
        );
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::kappa() const
{
    return alpha1()*thermo1_->kappa() + alpha2()*thermo2_->kappa();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::kappa
(
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappa(patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappa(patchi);
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::alphahe() const
{
    return
        alpha1()*thermo1_->alphahe()
      + alpha2()*thermo2_->alphahe();
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::alphahe
(
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->alphahe(patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->alphahe(patchi);
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->kappaEff(alphat)
      + alpha2()*thermo2_->kappaEff(alphat);
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappaEff(alphat, patchi);
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->alphaEff(alphat)
      + alpha2()*thermo2_->alphaEff(alphat);
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->alphaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->alphaEff(alphat, patchi);
}
bool Foam::twoPhaseMixtureThermo::read()
{
    if (psiThermo::read())
    {
        return interfaceProperties::read();
    }

    return false;
}
// ************************************************************************* //
