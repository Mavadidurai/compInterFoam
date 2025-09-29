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
#include "autoPtr.H"

#include "Pstream.H"
#include "Switch.H"
#include "Tuple2.H"
#include "Time.H"
#include <cmath>
#include <string>
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam { defineTypeNameAndDebug(twoPhaseMixtureThermo, 0); }
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
    thermo2_(nullptr),  // Initialise phase-change properties; values populated from metal dictionary
    latentHeat_(0.0),
    T_melt_(0.0),
    T_vapor_(0.0),
    dtFloor_(1e-12),
    Q_laser_
    (
        IOobject("Q_laser", U.mesh().time().timeName(), U.mesh(), IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        U.mesh(),
        dimensionedScalar("Q0", dimPower/dimVolume, 0.0)  // [W/m³]
    ),
    phaseChangeSource_
    (
        IOobject("phaseChangeSource", U.mesh().time().timeName(), U.mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
        U.mesh(),
        dimensionedScalar("source", dimTemperature/dimTime, 0.0)  // [K/s]
    ),
    phaseChangeRelaxCoeff_
    (
        IOobject("phaseChangeRelaxCoeff", U.mesh().time().timeName(), U.mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
        U.mesh(),
        dimensionedScalar("relax", dimless/dimTime, 0.0)
    ),
    dgdt_
    (
        IOobject("dgdt", U.mesh().time().timeName(), U.mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
        U.mesh(),
        dimensionedScalar("dgdt", dimless/dimTime, 0.0)
    ),
    nu1_("nu1", dimViscosity, 0.0),
    nu2_("nu2", dimViscosity, 0.0),
    rho1_("rho1", dimDensity, 0.0),
    rho2_("rho2", dimDensity, 0.0),
    ClTTM_("ClTTM", dimEnergy/dimVolume/dimTemperature, 0.0)
{
    const dictionary& transportDict =
        U.mesh().lookupObject<dictionary>("transportProperties");
    const dictionary& controlDict = U.mesh().time().controlDict();

    auto tryLookup =
        [&](const dictionary& dict, const word& entryName, const char* location, scalar& value) -> bool
        {
            if (!dict.found(entryName))
            {
                return false;
            }

            dimensionedScalar dimValue(entryName, dict);
            value = dimValue.value();

            if (!std::isfinite(value))
            {
                FatalIOErrorInFunction(dict)
                    << "Entry '" << entryName
                    << "' in " << location << " is not finite"
                    << exit(FatalIOError);
            }

            return true;
        };

    const dictionary* phaseChangeDictPtr = nullptr;
    const char* phaseChangeDictLocation = nullptr;
    if (controlDict.found("phaseChangeCoeffs"))
    {
        phaseChangeDictPtr = &controlDict.subDict("phaseChangeCoeffs");
        phaseChangeDictLocation = "controlDict.phaseChangeCoeffs";
    }
    else
    {
        FatalIOErrorInFunction(controlDict)
            << "Missing required dictionary 'phaseChangeCoeffs' in controlDict"
            << " (controlDict.phaseChangeCoeffs)"
            << exit(FatalIOError);
    }

    auto lookupOptionalScalar =
        [&](const word& entryName, scalar& value) -> bool
        {
            auto tryNameInAllLocations =
                [&](const word& name) -> bool
                {
                    if
                    (
                        phaseChangeDictPtr
                     && tryLookup
                        (
                            *phaseChangeDictPtr,
                            name,
                            phaseChangeDictLocation
                                ? phaseChangeDictLocation
                                : "phaseChangeCoeffs",
                            value
                        )
                    )
                    {
                        return true;
                    }

                    if
                    (
                        tryLookup
                        (
                            controlDict,
                            name,
                            "controlDict",
                            value
                        )
                    )
                    {
                        return true;
                    }

                    return false;
                };

            if (tryNameInAllLocations(entryName))
            {
                return true;
            }

            if
            (
                entryName == "T_vapor"
             || entryName == "Tvapor"
             || entryName == "Tvap"
            )
            {
                if (entryName != "Tvapor" && tryNameInAllLocations(word("Tvapor")))
                {
                    return true;
                }
                if (entryName != "T_vapor" && tryNameInAllLocations(word("T_vapor")))
                {
                    return true;
                }
                if (entryName != "Tvap" && tryNameInAllLocations(word("Tvap")))
                {
                    return true;
                }
            }

            return false;
        };

    auto lookupRequiredScalar =
        [&](const word& primaryName, const word& legacyName) -> scalar
        {
            scalar value = 0.0;

            if (!primaryName.empty() && lookupOptionalScalar(primaryName, value))
            {
                return value;
            }

            if (!legacyName.empty() && lookupOptionalScalar(legacyName, value))
            {
                return value;
            }

            FatalIOErrorInFunction(controlDict)
                << "Missing required entry '" << primaryName
                << "' (or legacy '" << legacyName
                << "') in controlDict or controlDict.phaseChangeCoeffs"
                << exit(FatalIOError);

            return 0.0;
        };

    latentHeat_ = lookupRequiredScalar("latentHeat", "hf");
    T_melt_ = lookupRequiredScalar("T_melt", "Tsol");
    T_vapor_ = lookupRequiredScalar("T_vapor", "Tvap");

    if (latentHeat_ <= SMALL)
    {
        FatalIOErrorInFunction(controlDict)
            << "Latent heat ('latentHeat' or legacy 'hf') in controlDict or"
            << " controlDict.phaseChangeCoeffs must be positive"
            << exit(FatalIOError);
    }

    if (T_melt_ <= 0 || T_vapor_ <= 0)
    {
        FatalIOErrorInFunction(controlDict)
            << "Phase change temperatures ('T_melt'/'T_vapor' or legacy 'Tsol'/'Tvap')"
            << " in controlDict or controlDict.phaseChangeCoeffs must be positive"
            << exit(FatalIOError);
    }

    if (T_melt_ >= T_vapor_)
    {
        FatalIOErrorInFunction(controlDict)
            << "Expected T_melt < T_vapor (legacy Tsol < Tvap) in controlDict"
            << " or controlDict.phaseChangeCoeffs"
            << exit(FatalIOError);
    }
    if (debug)
    {
        const Switch writePhaseTemperatures =
            U.mesh().time().controlDict().lookupOrDefault<Switch>
            (
                "writePhaseTemperatures",
                false
            );

        if (writePhaseTemperatures)
        {
            T_.write();
            fileHandler().flush();
        }
    }

    auto readDimensionedScalar =
        [&](const dictionary& dict,
            const word& entryName,
            const std::string& location) -> dimensionedScalar
        {
            if (!dict.found(entryName))
            {
                FatalIOErrorInFunction(dict)
                    << "Missing required entry '" << entryName
                    << "' in " << location
                    << exit(FatalIOError);
            }

            dimensionedScalar value(entryName, dict);
            const scalar val = value.value();

            if (!std::isfinite(val) || val <= 0)
            {
                FatalIOErrorInFunction(dict)
                    << "Entry '" << entryName << "' in " << location
                    << " must be finite and positive"
                    << exit(FatalIOError);
            }

            return value;
        };

    if
    (
        transportDict.found(phase1Name())
     && transportDict.found(phase2Name())
     && transportDict.isDict(phase1Name())
     && transportDict.isDict(phase2Name())
    )
    {
        const dictionary& phase1Dict = transportDict.subDict(phase1Name());
        const dictionary& phase2Dict = transportDict.subDict(phase2Name());

        nu1_ = readDimensionedScalar
        (
            phase1Dict,
            "nu",
            std::string("transportProperties.") + phase1Name()
        );
        nu2_ = readDimensionedScalar
        (
            phase2Dict,
            "nu",
            std::string("transportProperties.") + phase2Name()
        );
        rho1_ = readDimensionedScalar
        (
            phase1Dict,
            "rho",
            std::string("transportProperties.") + phase1Name()
        );
        rho2_ = readDimensionedScalar
        (
            phase2Dict,
            "rho",
            std::string("transportProperties.") + phase2Name()
        );
    }
    else
    {
        nu1_ = readDimensionedScalar
        (
            transportDict,
            "nu1",
            "transportProperties"
        );
        nu2_ = readDimensionedScalar
        (
            transportDict,
            "nu2",
            "transportProperties"
        );
        rho1_ = readDimensionedScalar
        (
            transportDict,
            "rho1",
            "transportProperties"
        );
        rho2_ = readDimensionedScalar
        (
            transportDict,
            "rho2",
            "transportProperties"
        );
    }
    if (Pstream::master())
    {
        Info<< "Transport properties:" << nl
            << "    " << phase1Name() << ": nu=" << nu1_.value()
            << ", rho=" << rho1_.value() << nl
            << "    " << phase2Name() << ": nu=" << nu2_.value()
            << ", rho=" << rho2_.value() << endl;
    }
    updateTwoTemperatureCache();
    thermo1_ = rhoThermo::New(U.mesh());
    thermo2_ = rhoThermo::New(U.mesh());
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
    computePhaseChange();
    phaseChangeSource_.correctBoundaryConditions();
    phaseChangeRelaxCoeff_.correctBoundaryConditions();
    // Compute mass-transfer rate [1/s]
    dgdt_ = computeMassTransfer()();
    dgdt_.correctBoundaryConditions();
   if (debug)
    {
        const fvMesh& mesh = T_.mesh();
        const label nCells = mesh.nCells();
        const label nSamples = Foam::min(label(5), nCells);
        if (nSamples > 0)
        {
            labelList sampleCells(nSamples);
            for (label i = 0; i < nSamples; ++i)
            {
                sampleCells[i] = i;
            }
            scalarField pSamples(nSamples, 0.0);
            scalarField tSamples(nSamples, 0.0);
            const scalarField& pInternal = p_.internalField();
            const scalarField& tInternal = T_.internalField();
            forAll(sampleCells, i)
            {
                const label celli = sampleCells[i];
                pSamples[i] = pInternal[celli];
                tSamples[i] = tInternal[celli];
            }
            tmp<scalarField> hSamplesTmp = he(pSamples, tSamples, sampleCells);
            const scalarField& hSamples = hSamplesTmp();
            tmp<scalarField> recoveredTmp = THE(hSamples, pSamples, tSamples, sampleCells);
            const scalarField& recovered = recoveredTmp();
            scalar maxError = 0.0;
            forAll(recovered, i)
            {
                maxError = Foam::max
                (
                    maxError,
                    Foam::mag(recovered[i] - tSamples[i])
                );
            }
            const scalar tolerance = 1e-6;
            if (maxError > tolerance)
            {
                FatalErrorInFunction
                    << "Regression check failed: THE(he(...)) deviates from"
                    << " the supplied temperature by " << maxError << " K"
                    << nl << "    tolerance = " << tolerance << " K" << nl
                    << exit(FatalError);
            }
        }
    }
}
void Foam::twoPhaseMixtureThermo::setQLaser(const volScalarField& src)
{
    Q_laser_ = src;
    Q_laser_.correctBoundaryConditions();
}
void Foam::twoPhaseMixtureThermo::updateTwoTemperatureCache()
{
    const Time& time = T_.mesh().time();
    const bool master = Pstream::master();
    if (time.controlDict().found("twoTemperatureProperties"))
    {
        const dictionary& twoTempDict =
            time.controlDict().subDict("twoTemperatureProperties");
        ClTTM_ = twoTempDict.lookupOrDefault<dimensionedScalar>("Cl", ClTTM_);
    }
    else
    {
        static bool warnedMissingTwoTemp = false;
        const scalar fallbackCl = 2.5e6;
        if (!warnedMissingTwoTemp)
        {
            if (master)
            {
                WarningInFunction
                    << "twoTemperatureProperties not found in controlDict; using"
                    << " fallback lattice heat capacity ClTTM=" << fallbackCl
                    << " J/m^3/K." << endl;
            }
            warnedMissingTwoTemp = true;
        }
        ClTTM_ = dimensionedScalar
        (
            "ClTTM",
            dimEnergy/dimVolume/dimTemperature,
            fallbackCl
        );
    }
}
void Foam::twoPhaseMixtureThermo::setClTTM(const dimensionedScalar& Cl)
{
    ClTTM_ = Cl;
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
        dimEnergy/dimMass,
        latentHeat_
    );
}
Foam::tmp<Foam::volScalarField>
Foam::twoPhaseMixtureThermo::sigma() const
{
    Foam::tmp<Foam::volScalarField> tsigma
    (
        surfaceTensionModel::New(*this, alpha1().mesh())->sigma()
    );
    if
    (
        tsigma.valid()
     && tsigma().dimensions() != dimForce/dimLength
    )
    {
        FatalErrorInFunction
            << "Surface tension field has dimensions "
            << tsigma().dimensions()
            << ", expected " << dimForce/dimLength << nl
            << exit(FatalError);
    }
    return tsigma;
}
void Foam::twoPhaseMixtureThermo::computePhaseChange()
{
    phaseChangeSource_ = dimensionedScalar("source", dimTemperature/dimTime, 0.0);
    // Reset the implicit relaxation coefficient
    phaseChangeRelaxCoeff_ = dimensionedScalar("relax", dimless/dimTime, 0.0);
    const bool master = Pstream::master();
    const fvMesh& mesh = T_.mesh();
    // Access coefficients from transportProperties
    const dictionary& transportDict = mesh.lookupObject<dictionary>("transportProperties");
    const dictionary& controlDict = mesh.time().controlDict();

    const volScalarField* TlPtr = nullptr;
    autoPtr<volScalarField> TlTmp;

    if (mesh.foundObject<volScalarField>("Tl"))
    {
        TlPtr = &mesh.lookupObject<volScalarField>("Tl");
    }
    else
    {
        IOobject TlHeader
        (
            "Tl",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (TlHeader.typeHeaderOk<volScalarField>())
        {
            TlTmp.reset(new volScalarField(TlHeader, mesh));
            TlPtr = TlTmp.ptr();
        }
    }

    if (!TlPtr)
    {
        static bool warnedMissingTl = false;
        if (!warnedMissingTl && master)
        {
            WarningInFunction
                << "Lattice temperature field 'Tl' not available;"
                << " skipping phase-change update until it is created." << endl;
            warnedMissingTl = true;
        }

        return;
    }

    const volScalarField& Tl = *TlPtr;
    const dictionary* pcPtr = nullptr;
    word pcLocation("phaseChangeCoeffs");
    if (controlDict.found("phaseChangeCoeffs"))
    {
        pcPtr = &controlDict.subDict("phaseChangeCoeffs");
        pcLocation = "controlDict.phaseChangeCoeffs";
    }
    else if (transportDict.found("phaseChangeCoeffs"))
    {
        pcPtr = &transportDict.subDict("phaseChangeCoeffs");
        pcLocation = "transportProperties.phaseChangeCoeffs";
    }

    if (!pcPtr)
    {
        FatalErrorInFunction
            << "phaseChangeCoeffs not found in controlDict or transportProperties" << nl
            << "Supply a phaseChangeCoeffs sub-dictionary or disable phase change." << nl
            << exit(FatalError);

    }
    const dictionary& pc = *pcPtr;
    const scalar Tvapor = pc.lookupOrDefault<scalar>("Tvapor", T_vapor_);
    const scalar windowWidth = pc.lookupOrDefault<scalar>("windowWidth", 0.0);
    dtFloor_ = pc.lookupOrDefault<scalar>("dtFloor", dtFloor_);
    scalar relaxationRate = pc.lookupOrDefault<scalar>("relaxationRate", -1.0);
    if (relaxationRate < 0)
    {
        const scalar relaxationTime = pc.lookupOrDefault<scalar>("relaxationTime", -1.0);
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
    const word maxSourceDefaultKey("phaseChangeMaxSourceDefault");
    const bool hasTransportDefault = transportDict.found(maxSourceDefaultKey);
    scalar defaultMaxSource = transportDict.lookupOrDefault<scalar>(maxSourceDefaultKey, 1e7);
    scalar maxSource = defaultMaxSource;
    bool usingDefaultMaxSource = false;
    if (pc.found("maxSource"))
    {
        maxSource = pc.lookupOrDefault<scalar>("maxSource", defaultMaxSource);
    }
    else if (pc.found("minCoefficient"))
    {
        maxSource = readScalar(pc.lookup("minCoefficient"));
    }
    else
    {
        usingDefaultMaxSource = true;
    }
    if (usingDefaultMaxSource)
    {
        static bool defaultMaxSourceWarned = false;
        if (!defaultMaxSourceWarned)
        {
            if (master)
            {
                WarningInFunction
                    << pcLocation << ":maxSource not specified; using "
                    << (hasTransportDefault
                        ? "transportProperties entry 'phaseChangeMaxSourceDefault'"
                        : "internal fallback")
                    << " = " << maxSource << " [K/s]" << endl;
            }
            defaultMaxSourceWarned = true;
        }
    }
    const bool limitSource = (maxSource > 0 && maxSource < GREAT);
    scalar metalCutoff = -1.0;
    if (pc.found("phaseChangeMetalCutoff"))
    {
        metalCutoff = pc.lookupOrDefault<scalar>("phaseChangeMetalCutoff", metalCutoff);
    }
    else if (transportDict.found("phaseChangeMetalCutoff"))
    {
        metalCutoff = transportDict.lookupOrDefault<scalar>
        (
            "phaseChangeMetalCutoff",
            metalCutoff
        );
    }
    if (metalCutoff < 0)
    {
        const Time& time = T_.mesh().time();
        if (time.controlDict().found("twoTemperatureProperties"))
        {
            const dictionary& twoTempDict =
                time.controlDict().subDict("twoTemperatureProperties");
            const scalar metalFloor =
                twoTempDict.lookupOrDefault<scalar>("metalFractionFloor", 1e-6);
            metalCutoff = twoTempDict.lookupOrDefault<scalar>
            (
                "metalFractionCutoff",
                metalFloor
            );
        }
        else
        {
            metalCutoff = 1e-6;
        }
    }
    metalCutoff = Foam::max(metalCutoff, 0.0);    
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
        if (master)
        {
            Info<< "phaseChangeCoeffs:" << nl
                << "    Tvapor        " << Tvapor << nl
                << "    windowWidth   " << windowWidth << nl
                << "    dtFloor       " << dtFloor_ << nl
                << "    relaxationRate " << relaxationRate << nl
                << "    relaxationTime "
                << (relaxationRate > SMALL ? 1.0/relaxationRate : GREAT) << nl
                << "    maxSource     " << maxSource << nl
                << "    metalCutoff   " << metalCutoff << nl
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
            if (master)
            {
                Info<< "phaseChangeCoeffs inactive at time " << timeVal << nl;
            }
            return;
        }
    }
    const scalar minAlpha = 0.0;
    const scalar maxAlpha = 1.0;
    const scalar localRateMax = 1.0/Foam::max(dtFloor_, SMALL);    
    forAll(Tl, cellI)
    {
        const scalar a1 = Foam::min(Foam::max(alpha1()[cellI], minAlpha), maxAlpha);
        const scalar Tcell = Tl[cellI];
        const bool metalActive = (a1 > metalCutoff);
        scalar available = 0.0;
        if (metalActive)
        {
            if (Tcell > Tvapor)
            {
                available = a1;
            }
            else if (Tcell < Tvapor)
            {
                available = 1.0 - a1;
            }
        }
        scalar localRate = relaxationRate*available;
        if (windowWidth > SMALL)
        {
            localRate *= Foam::min(Foam::mag(Tcell - Tvapor)/windowWidth, 1.0);
        }
        scalar sourceVal = localRate*(Tvapor - Tcell);
        if (limitSource)
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
        if (!metalActive)
        {
            localRate = 0.0;
            sourceVal = 0.0;
        }
        localRate = Foam::min(localRate, localRateMax);
        phaseChangeRelaxCoeff_[cellI] = localRate;
        phaseChangeSource_[cellI] = sourceVal;
    }
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
                IOobject::NO_WRITE,
                false
            ),
            T_.mesh(),
            dimensionedScalar("dgdt", dimless/dimTime, 0.0)
        )
    );
    volScalarField& dgdt = tDgdt.ref();
    const dictionary& transportDict =
        T_.mesh().lookupObject<dictionary>("transportProperties");
    const dictionary& controlDict = T_.mesh().time().controlDict();
    const bool master = Pstream::master();
    const dictionary* mtPtr = nullptr;
    word mtLocation("massTransferCoeffs");
    if (controlDict.found("massTransferCoeffs"))
    {
        mtPtr = &controlDict.subDict("massTransferCoeffs");
        mtLocation = "controlDict.massTransferCoeffs";
    }
    else if (transportDict.found("massTransferCoeffs"))
    {
        mtPtr = &transportDict.subDict("massTransferCoeffs");
        mtLocation = "transportProperties.massTransferCoeffs";
    }

    if (!mtPtr)
    {
        if (master)
        {
            Info<< "massTransferCoeffs not found - skipping mass transfer" << nl;
        }
        return tDgdt;
    }
    const dictionary& mt = *mtPtr;
    const scalar rateMax = mt.lookupOrDefault<scalar>("rateMax", -1.0);
    // lattice heat capacity [J/m^3/K] from two-temperature properties
    const scalar ClVal = ClTTM_.value();
    const scalar LVal = latentHeat().value();
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
        if (master)
        {
            Info<< mtLocation << ':' << nl;
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
        }
        loggedMassTransfer = true;
    }
    if (ClVal <= SMALL)
    {
        if (master)
        {
            WarningInFunction
                << "Cached lattice heat capacity ClTTM=" << ClVal
                << " J/m^3/K is non-positive; skipping mass transfer computation."
                << endl;
        }
        return tDgdt;
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
            if (master)
            {
                Info<< "massTransferCoeffs inactive at time " << timeVal << nl;
            }
            return tDgdt;
        }
    }
    // Retrieve phase-1 density field and store temporary to avoid referencing
    // destroyed objects
    const tmp<volScalarField> rho1Tmp = thermo1_->rho();
    const volScalarField& rho1Field = rho1Tmp();
    scalar rho1Ref = rho1_.value();
    if (rho1Ref <= SMALL)
    {
        rho1Ref = Foam::gMax(rho1Field.internalField());
    }
    const scalar rho1Min = Foam::max(1e-6*rho1Ref, SMALL);
    forAll(dgdt, cellI)
    {
        const scalar source = phaseChangeSource_[cellI];
        if (Foam::mag(source) <= SMALL)
        {
            dgdt[cellI] = 0.0;
            continue;
        }
        scalar localRate =
            -(ClVal*source)
            /(LVal*max(rho1Field[cellI], rho1Min));
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
namespace
{
    inline Foam::scalar massWeighted
    (
        const Foam::scalar alpha1,
        const Foam::scalar rho1,
        const Foam::scalar value1,
        const Foam::scalar alpha2,
        const Foam::scalar rho2,
        const Foam::scalar value2
    )
    {
        const Foam::scalar rhoMix =
            Foam::max(alpha1*rho1 + alpha2*rho2, Foam::SMALL);
        return
            (alpha1*rho1*value1 + alpha2*rho2*value2)
           /rhoMix;
    }
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const volScalarField& alpha1Field = alpha1();
    const volScalarField& alpha2Field = alpha2();
    tmp<volScalarField> the1 = thermo1_->he(p, T);
    tmp<volScalarField> the2 = thermo2_->he(p, T);
    tmp<volScalarField> trho1 = thermo1_->rho();
    tmp<volScalarField> trho2 = thermo2_->rho();
    tmp<volScalarField> numerator =
        alpha1Field*trho1()*the1()
      + alpha2Field*trho2()*the2();
    tmp<volScalarField> denominator =
        alpha1Field*trho1()
      + alpha2Field*trho2();
    return numerator
       /
        (
            denominator
          + dimensionedScalar("rhoMin", dimDensity, Foam::SMALL)
        );
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    const label nCells = cells.size();
    scalarField alpha1Field(nCells, 0.0);
    scalarField alpha2Field(nCells, 0.0);
    forAll(cells, i)
    {
        const label celli = cells[i];
        alpha1Field[i] = alpha1()[celli];
        alpha2Field[i] = alpha2()[celli];
    }
    tmp<scalarField> the1 = thermo1_->he(p, T, cells);
    tmp<scalarField> the2 = thermo2_->he(p, T, cells);
    tmp<scalarField> trho1 = thermo1_->rhoEoS(p, T, cells);
    tmp<scalarField> trho2 = thermo2_->rhoEoS(p, T, cells);
    scalarField result(nCells, 0.0);
    forAll(result, i)
    {
        result[i] = massWeighted
        (
            alpha1Field[i],
            trho1()[i],
            the1()[i],
            alpha2Field[i],
            trho2()[i],
            the2()[i]
        );
    }
    return tmp<scalarField>(new scalarField(std::move(result)));
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const fvPatchScalarField& alpha1Patch = alpha1().boundaryField()[patchi];
    const fvPatchScalarField& alpha2Patch = alpha2().boundaryField()[patchi];
    tmp<scalarField> the1 = thermo1_->he(p, T, patchi);
    tmp<scalarField> the2 = thermo2_->he(p, T, patchi);
    tmp<scalarField> trho1 = thermo1_->rho(patchi);
    tmp<scalarField> trho2 = thermo2_->rho(patchi);
    scalarField result(the1().size(), 0.0);
    forAll(result, facei)
    {
        result[facei] = massWeighted
        (
            alpha1Patch[facei],
            trho1()[facei],
            the1()[facei],
            alpha2Patch[facei],
            trho2()[facei],
            the2()[facei]
        );
    }
    return tmp<scalarField>(new scalarField(std::move(result)));
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hc() const
{
    return alpha1()*thermo1_->hc() + alpha2()*thermo2_->hc();
}
namespace
{
    template<class RhoEval1, class RhoEval2, class HeEval1, class HeEval2, class CpEval1, class CpEval2>
    Foam::scalarField invertMassWeightedEnthalpy
    (
        const Foam::scalarField& h,
        const Foam::scalarField& p,
        const Foam::scalarField& T0,
        const Foam::scalarField& alpha1,
        const Foam::scalarField& alpha2,
        const RhoEval1& rhoEval1,
        const RhoEval2& rhoEval2,
        const HeEval1& heEval1,
        const HeEval2& heEval2,
        const CpEval1& cpEval1,
        const CpEval2& cpEval2
    )
    {
        Foam::scalarField T(T0);
        const Foam::label maxIter = 50;
        const Foam::scalar tol = 1e-6;
        bool converged = false;
        for (Foam::label iter = 0; iter < maxIter; ++iter)
        {
            Foam::tmp<Foam::scalarField> trho1 = rhoEval1(p, T);
            Foam::tmp<Foam::scalarField> trho2 = rhoEval2(p, T);
            Foam::tmp<Foam::scalarField> the1 = heEval1(p, T);
            Foam::tmp<Foam::scalarField> the2 = heEval2(p, T);
            Foam::tmp<Foam::scalarField> tcp1 = cpEval1(p, T);
            Foam::tmp<Foam::scalarField> tcp2 = cpEval2(p, T);
            Foam::scalar maxDeltaT = 0.0;
            forAll(T, i)
            {
                const Foam::scalar rho1 = trho1()[i];
                const Foam::scalar rho2 = trho2()[i];
                const Foam::scalar he1 = the1()[i];
                const Foam::scalar he2 = the2()[i];
                const Foam::scalar cp1 = tcp1()[i];
                const Foam::scalar cp2 = tcp2()[i];
                const Foam::scalar a1 = alpha1[i];
                const Foam::scalar a2 = alpha2[i];
                const Foam::scalar rhoMix = Foam::max(a1*rho1 + a2*rho2, Foam::SMALL);
                const Foam::scalar mixtureHe =
                    (a1*rho1*he1 + a2*rho2*he2)/rhoMix;
                const Foam::scalar mixtureCp =
                    (a1*rho1*cp1 + a2*rho2*cp2)/rhoMix;
                const Foam::scalar deltaT =
                    (mixtureHe - h[i])
                   /Foam::max(mixtureCp, Foam::SMALL);
                T[i] -= deltaT;
                maxDeltaT = Foam::max(maxDeltaT, Foam::mag(deltaT));
            }
            if (maxDeltaT < tol)
            {
                converged = true;
                break;
            }
        }
        if (!converged)
        {
            WarningInFunction
                << "THE() iteration reached the maximum number of iterations"
                << " without achieving the requested tolerance." << Foam::endl;
        }
        return T;
    }
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    const label nCells = cells.size();
    scalarField alpha1Field(nCells, 0.0);
    scalarField alpha2Field(nCells, 0.0);
    forAll(cells, i)
    {
        const label celli = cells[i];
        alpha1Field[i] = alpha1()[celli];
        alpha2Field[i] = alpha2()[celli];
    }
    auto rhoEval1 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo1_->rhoEoS(pVals, TVals, cells);
    };
    auto rhoEval2 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo2_->rhoEoS(pVals, TVals, cells);
    };
    auto heEval1 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo1_->he(pVals, TVals, cells);
    };
    auto heEval2 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo2_->he(pVals, TVals, cells);
    };
    auto cpEval1 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo1_->Cp(pVals, TVals, cells);
    };
    auto cpEval2 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo2_->Cp(pVals, TVals, cells);
    };
    scalarField result = invertMassWeightedEnthalpy
    (
        h,
        p,
        T0,
        alpha1Field,
        alpha2Field,
        rhoEval1,
        rhoEval2,
        heEval1,
        heEval2,
        cpEval1,
        cpEval2
    );
    return tmp<scalarField>(new scalarField(std::move(result)));
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    const fvPatchScalarField& alpha1Patch = alpha1().boundaryField()[patchi];
    const fvPatchScalarField& alpha2Patch = alpha2().boundaryField()[patchi];
    scalarField alpha1Field(alpha1Patch.size(), 0.0);
    scalarField alpha2Field(alpha2Patch.size(), 0.0);
    forAll(alpha1Field, facei)
    {
        alpha1Field[facei] = alpha1Patch[facei];
        alpha2Field[facei] = alpha2Patch[facei];
    }
    auto rhoEval1 = [&](const scalarField&, const scalarField&)
    {
        return thermo1_->rho(patchi);
    };
    auto rhoEval2 = [&](const scalarField&, const scalarField&)
    {
        return thermo2_->rho(patchi);
    };
    auto heEval1 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo1_->he(pVals, TVals, patchi);
    };
    auto heEval2 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo2_->he(pVals, TVals, patchi);
    };
    auto cpEval1 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo1_->Cp(pVals, TVals, patchi);
    };
    auto cpEval2 = [&](const scalarField& pVals, const scalarField& TVals)
    {
        return thermo2_->Cp(pVals, TVals, patchi);
    };
    scalarField result = invertMassWeightedEnthalpy
    (
        h,
        p,
        T0,
        alpha1Field,
        alpha2Field,
        rhoEval1,
        rhoEval2,
        heEval1,
        heEval2,
        cpEval1,
        cpEval2
    );
    return tmp<scalarField>(new scalarField(std::move(result)));
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
        updateTwoTemperatureCache();        
        return interfaceProperties::read();
    }
    return false;
}
// ************************************************************************* //
