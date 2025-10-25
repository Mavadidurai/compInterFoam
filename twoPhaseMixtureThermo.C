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
    gasConstant_(8.314/0.04787),
    evaporationCoeff_(0.18),
    relaxationTime_(1e-12),
    alphaMin_(0.01),
    alphaMax_(1.0),
    onlyAboveVapor_(false),
    activationWindows_(),
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
    phaseChangeMassFlux_
    (
        IOobject
        (
            "phaseChangeMassFlux",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar
        (
            "massFlux",
            dimensionSet(1, -2, -1, 0, 0, 0, 0),
            0.0
        )
    ),
    nu1_("nu1", dimViscosity, 0.0),
    nu2_("nu2", dimViscosity, 0.0),
    rho1_("rho1", dimDensity, 0.0),
    rho2_("rho2", dimDensity, 0.0),
    ClTTM_("ClTTM", dimEnergy/dimVolume/dimTemperature, 0.0),
    sigmaModel_(nullptr)
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
    std::string phaseChangeDictLocation;
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

    const std::string phaseChangeDictDisplay =
        phaseChangeDictLocation.empty()
            ? "controlDict.phaseChangeCoeffs"
            : phaseChangeDictLocation;
    const std::string controlDictDisplay = "controlDict";

    auto isPhaseChangeCoeffsLocation =
        [&](const std::string& location) -> bool
        {
            return !location.empty()
                && location.find(phaseChangeDictDisplay) == 0;
        };

    auto isControlDictLocation =
        [&](const std::string& location) -> bool
        {
            return !location.empty()
                && location.find(controlDictDisplay) == 0
                && !isPhaseChangeCoeffsLocation(location);
        };

    auto alternateDictionary =
        [&](const std::string& location) -> const std::string&
        {
            return isControlDictLocation(location)
                ? phaseChangeDictDisplay
                : controlDictDisplay;
        };

    auto combinedAlternateDictionary =
        [&](const std::string& a, const std::string& b) -> std::string
        {
            const bool aIsPhase = isPhaseChangeCoeffsLocation(a);
            const bool bIsPhase = isPhaseChangeCoeffsLocation(b);
            const bool aIsControl = isControlDictLocation(a);
            const bool bIsControl = isControlDictLocation(b);

            if ((aIsPhase || bIsPhase) && (aIsControl || bIsControl))
            {
                return phaseChangeDictDisplay + " and " + controlDictDisplay;
            }

            if (aIsPhase || bIsPhase)
            {
                return controlDictDisplay;
            }

            return phaseChangeDictDisplay;
        };

    auto lookupOptionalScalar =
        [&](const word& entryName, scalar& value, std::string& location) -> bool
        {
            auto tryNameInAllLocations =
                [&](const word& name, const bool legacyAlias) -> bool
                {
                    if
                    (
                        phaseChangeDictPtr
                     && tryLookup
                        (
                            *phaseChangeDictPtr,
                            name,
                            phaseChangeDictDisplay.c_str(),
                            value
                        )
                    )
                    {
                        location = phaseChangeDictDisplay;
                        if (legacyAlias)
                        {
                            location += " (legacy entry '" + name + "')";
                        }
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
                        location = "controlDict";
                        if (legacyAlias)
                        {
                            location += " (legacy entry '" + name + "')";
                        }
                        return true;
                    }

                    return false;
                };

            if (tryNameInAllLocations(entryName, false))
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
                if
                (
                    entryName != "Tvapor"
                 && tryNameInAllLocations(word("Tvapor"), true)
                )
                {
                    return true;
                }
                if
                (
                    entryName != "T_vapor"
                 && tryNameInAllLocations(word("T_vapor"), true)
                )
                {
                    return true;
                }
                if
                (
                    entryName != "Tvap"
                 && tryNameInAllLocations(word("Tvap"), true)
                )
                {
                    return true;
                }
            }

            return false;
        };

    auto lookupRequiredScalar =
        [&](const word& primaryName, const word& legacyName, std::string& location) -> scalar
        {
            scalar value = 0.0;

            if
            (
                !primaryName.empty()
             && lookupOptionalScalar(primaryName, value, location)
            )
            {
                return value;
            }

            if
            (
                !legacyName.empty()
             && lookupOptionalScalar(legacyName, value, location)
            )
            {
                location += " (legacy entry '" + legacyName + "')";
                return value;
            }

            FatalIOErrorInFunction(controlDict)
                << "Missing required entry '" << primaryName
                << "' (or legacy '" << legacyName
                << "') in " << phaseChangeDictDisplay << " or controlDict"
                << exit(FatalIOError);

            return 0.0;
        };

    std::string latentHeatLocation;
    std::string TmeltLocation;
    std::string TvaporLocation;

    latentHeat_ = lookupRequiredScalar("latentHeat", "hf", latentHeatLocation);
    T_melt_ = lookupRequiredScalar("T_melt", "Tsol", TmeltLocation);
    T_vapor_ = lookupRequiredScalar("T_vapor", "Tvap", TvaporLocation);

    if (latentHeat_ <= SMALL)
    {
        FatalIOErrorInFunction(controlDict)
            << "Phase-change latent heat ('latentHeat' or legacy 'hf') in "
            << latentHeatLocation
            << " must be positive (also checked "
            << alternateDictionary(latentHeatLocation)
            << ")"
            << exit(FatalIOError);
    }

    if (T_melt_ <= 0 || T_vapor_ <= 0)
    {
        FatalIOErrorInFunction(controlDict)
            << "Phase change temperatures ('T_melt'/'T_vapor' or legacy 'Tsol'/'Tvap')"
            << " in " << TmeltLocation << " and " << TvaporLocation
            << " must be positive (also checked "
            << combinedAlternateDictionary(TmeltLocation, TvaporLocation)
            << ")"
            << exit(FatalIOError);
    }

    if (T_melt_ >= T_vapor_)
    {
        FatalIOErrorInFunction(controlDict)
            << "Expected T_melt < T_vapor (legacy Tsol < Tvap) with values from "
            << TmeltLocation << " and " << TvaporLocation
            << " (also checked "
            << combinedAlternateDictionary(TmeltLocation, TvaporLocation)
            << ")"
            << exit(FatalIOError);
    }
    
    const dictionary& phaseChangeDict = *phaseChangeDictPtr;

    onlyAboveVapor_ = phaseChangeDict.lookupOrDefault<Switch>
    (
        "onlyAboveVapor",
        false
    );
    
    auto checkDimensions =
        [&](const word& entryName,
            const dimensionedScalar& value,
            const dimensionSet& expected) -> void
        {
            if (value.dimensions() != expected)
            {
                FatalIOErrorInFunction(phaseChangeDict)
                    << "Entry '" << entryName << "' in "
                    << phaseChangeDictDisplay
                    << " has dimensions " << value.dimensions()
                    << " but expected " << expected
                    << exit(FatalIOError);
            }
        };

    auto ensureFinite =
        [&](const word& entryName, const scalar value) -> void
        {
            if (!std::isfinite(value))
            {
                FatalIOErrorInFunction(phaseChangeDict)
                    << "Entry '" << entryName << "' in "
                    << phaseChangeDictDisplay
                    << " must be finite"
                    << exit(FatalIOError);
            }
        };
    if (phaseChangeDict.found("tStart") || phaseChangeDict.found("tEnd"))
    {
        if (!(phaseChangeDict.found("tStart") && phaseChangeDict.found("tEnd")))
        {
            FatalIOErrorInFunction(phaseChangeDict)
                << "Both 'tStart' and 'tEnd' must be provided in "
                << phaseChangeDictDisplay
                << " to configure activation windows"
                << exit(FatalIOError);
        }

        List<scalar> tStart;
        List<scalar> tEnd;
        phaseChangeDict.lookup("tStart") >> tStart;
        phaseChangeDict.lookup("tEnd") >> tEnd;

        if (tStart.size() != tEnd.size())
        {
            FatalIOErrorInFunction(phaseChangeDict)
                << "Entries 'tStart' and 'tEnd' in " << phaseChangeDictDisplay
                << " must have the same number of values"
                << exit(FatalIOError);
        }

        activationWindows_.setSize(tStart.size());
        scalar previousEnd = -GREAT;
        forAll(tStart, i)
        {
            const scalar start = tStart[i];
            const scalar end = tEnd[i];

            if (!std::isfinite(start) || !std::isfinite(end))
            {
                FatalIOErrorInFunction(phaseChangeDict)
                    << "Entries 'tStart'/'tEnd' in " << phaseChangeDictDisplay
                    << " must contain finite values"
                    << exit(FatalIOError);
            }

            if (end <= start)
            {
                FatalIOErrorInFunction(phaseChangeDict)
                    << "Activation window requires tEnd > tStart in "
                    << phaseChangeDictDisplay << " (window " << i << ')'
                    << exit(FatalIOError);
            }

            if (i > 0 && start < previousEnd)
            {
                FatalIOErrorInFunction(phaseChangeDict)
                    << "Activation windows in " << phaseChangeDictDisplay
                    << " must be ordered and non-overlapping"
                    << exit(FatalIOError);
            }

            activationWindows_[i] = Tuple2<scalar, scalar>(start, end);
            previousEnd = end;
        }
    }

    if (phaseChangeDict.found("gasConstant"))
    {
        const dimensionedScalar value("gasConstant", phaseChangeDict);
        checkDimensions("gasConstant", value, dimEnergy/dimMass/dimTemperature);
        const scalar val = value.value();
        ensureFinite("gasConstant", val);

        if (val <= SMALL)
        {
            FatalIOErrorInFunction(phaseChangeDict)
                << "Entry 'gasConstant' in " << phaseChangeDictDisplay
                << " must be positive"
                << exit(FatalIOError);
        }

        gasConstant_ = val;
    }

    if (phaseChangeDict.found("evaporationCoeff"))
    {
        const dimensionedScalar value("evaporationCoeff", phaseChangeDict);
        checkDimensions("evaporationCoeff", value, dimless);
        const scalar val = value.value();
        ensureFinite("evaporationCoeff", val);

        if (val <= SMALL)
        {
            FatalIOErrorInFunction(phaseChangeDict)
                << "Entry 'evaporationCoeff' in " << phaseChangeDictDisplay
                << " must be positive"
                << exit(FatalIOError);
        }

        evaporationCoeff_ = val;
    }

    if (phaseChangeDict.found("evapRelaxationTime"))
    {
        const dimensionedScalar value("evapRelaxationTime", phaseChangeDict);
        checkDimensions("evapRelaxationTime", value, dimTime);
        const scalar val = value.value();
        ensureFinite("evapRelaxationTime", val);

        if (val <= SMALL)
        {
            FatalIOErrorInFunction(phaseChangeDict)
                << "Entry 'evapRelaxationTime' in " << phaseChangeDictDisplay
                << " must be positive"
                << exit(FatalIOError);
        }

        relaxationTime_ = val;
    }

    if (phaseChangeDict.found("alphaMin"))
    {
        const dimensionedScalar value("alphaMin", phaseChangeDict);
        checkDimensions("alphaMin", value, dimless);
        const scalar val = value.value();
        ensureFinite("alphaMin", val);

        if (val < 0 || val >= 1)
        {
            FatalIOErrorInFunction(phaseChangeDict)
                << "Entry 'alphaMin' in " << phaseChangeDictDisplay
                << " must satisfy 0 <= alphaMin < 1"
                << exit(FatalIOError);
        }

        alphaMin_ = val;
    }

    if (phaseChangeDict.found("alphaMax"))
    {
        const dimensionedScalar value("alphaMax", phaseChangeDict);
        checkDimensions("alphaMax", value, dimless);
        const scalar val = value.value();
        ensureFinite("alphaMax", val);

        if (val <= 0 || val > 1)
        {
            FatalIOErrorInFunction(phaseChangeDict)
                << "Entry 'alphaMax' in " << phaseChangeDictDisplay
                << " must satisfy 0 < alphaMax <= 1"
                << exit(FatalIOError);
        }

        alphaMax_ = val;
    }

    if (alphaMin_ >= alphaMax_)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Expected alphaMin < alphaMax in " << phaseChangeDictDisplay
            << ". Received alphaMin=" << alphaMin_
            << " and alphaMax=" << alphaMax_
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
    try
    {
        thermo1_ = rhoThermo::New(U.mesh(), phase1Name());
    }
    catch (const Foam::IOerror&)
    {
        FatalErrorInFunction
            << "Failed to construct rhoThermo for phase '" << phase1Name()
            << "'. Expected dictionary 'thermophysicalProperties."
            << phase1Name() << "'." << exit(FatalError);
    }

    try
    {
        thermo2_ = rhoThermo::New(U.mesh(), phase2Name());
    }
    catch (const Foam::IOerror&)
    {
        FatalErrorInFunction
            << "Failed to construct rhoThermo for phase '" << phase2Name()
            << "'. Expected dictionary 'thermophysicalProperties."
            << phase2Name() << "'." << exit(FatalError);
    }
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
    sigmaModel_.reset(nullptr);
    psi_ = alpha1()*thermo1_->psi() + alpha2()*thermo2_->psi();
    mu_ = alpha1()*thermo1_->mu() + alpha2()*thermo2_->mu();
    alpha_ = alpha1()*thermo1_->alpha() + alpha2()*thermo2_->alpha();
    interfaceProperties::correct();
    // Laser heating is now supplied externally via setQLaser()
    computePhaseChange();
    phaseChangeSource_.correctBoundaryConditions();
    phaseChangeRelaxCoeff_.correctBoundaryConditions();
    phaseChangeMassFlux_.correctBoundaryConditions();
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
    if (!sigmaModel_.valid())
    {
        sigmaModel_ = surfaceTensionModel::New(*this, alpha1().mesh());
    }

    Foam::tmp<Foam::volScalarField> tsigma
    (
        sigmaModel_->sigma()
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
    phaseChangeSource_.primitiveFieldRef() = 0.0;
    phaseChangeRelaxCoeff_.primitiveFieldRef() = 0.0;
    phaseChangeMassFlux_.primitiveFieldRef() = 0.0;
    phaseChangeMassFlux_.boundaryFieldRef() = 0.0;

    const fvMesh& mesh = T_.mesh();
    const Time& time = mesh.time();
    const scalar currentTime = time.value();
    const scalar laserEnd =
        time.controlDict().lookupOrDefault<scalar>("laserEndTime", GREAT);
    const scalar phaseChangeGrace = 10.0*relaxationTime_;

    bool withinWindow = activationWindows_.empty();
    forAll(activationWindows_, windowI)
    {
        const Tuple2<scalar, scalar>& window = activationWindows_[windowI];
        if
        (
            currentTime >= window.first()
         && currentTime <= window.second()
        )
        {
            withinWindow = true;
            break;
        }
    }

    static bool loggedInactivePhaseChange = false;

    if (!withinWindow)
    {
        if (!loggedInactivePhaseChange && Pstream::master())
        {
            Info<< "Phase-change source inactive outside configured window at t="
                << currentTime << endl;
        }
        loggedInactivePhaseChange = true;
        return;
    }

    if (loggedInactivePhaseChange)
    {
        loggedInactivePhaseChange = false;
    }

    if (currentTime > laserEnd + phaseChangeGrace)
    {
        return;
    }
    
    const volScalarField* TlPtr = mesh.foundObject<volScalarField>("Tl")
        ? &mesh.lookupObject<volScalarField>("Tl")
        : nullptr;

    if (!TlPtr)
    {
        return;
    }

    const volScalarField& Tl = *TlPtr;
    const volScalarField& p = mesh.lookupObject<volScalarField>("p");

    const scalar L = latentHeat_;
    const scalar R = gasConstant_;
    const scalar T_vap = T_vapor_;
    const scalar p_ref = 101325;

    const scalarField& alpha1Field = alpha1().primitiveField();
    const scalarField& TlField = Tl.primitiveField();
    const scalarField& pField = p.primitiveField();

    const bool enforceUpper = alphaMax_ < (scalar(1) - Foam::SMALL);

    const scalarField& cellVolumes = mesh.V();
    const bool restrictToVapor = onlyAboveVapor_;

    forAll(TlField, cellI)
    {
        const scalar T_local = TlField[cellI];
        const scalar p_local = pField[cellI];
        const scalar alpha = alpha1Field[cellI];

        if (alpha < alphaMin_ || (enforceUpper && alpha > alphaMax_))
        {
            continue;
        }

        if (restrictToVapor && T_local < T_vap)
        {
            continue;
        }

        const scalar inv_T = 1.0/Foam::max(T_local, 100.0);
        const scalar inv_Tvap = 1.0/T_vap;
        const scalar p_vapor = p_ref*std::exp(L/R*(inv_Tvap - inv_T));

        const scalar sqrt_T = std::sqrt(Foam::max(T_local, 1.0));
        const scalar sqrt_2piR = std::sqrt(2.0*3.14159*R);

        const scalar j_evap = evaporationCoeff_*p_vapor/(sqrt_2piR*sqrt_T);
        const scalar j_cond = evaporationCoeff_*p_local/(sqrt_2piR*sqrt_T);

        const scalar j_net = j_evap - j_cond;

        const scalar Cl = ClTTM_.value();

        if (Cl <= SMALL)
        {
            continue;
        }

        const scalar heatFlux = j_net*L;  // [W/m^2]
        const scalar cellVolume = cellVolumes[cellI];

        if (cellVolume <= VSMALL)
        {
            continue;
        }

        const scalar metalFraction = Foam::max(alpha, scalar(0));

        if (metalFraction <= VSMALL)
        {
            continue;
        }

        const scalar cellLength = std::pow(cellVolume, 1.0/3.0);
        const scalar meltThickness = Foam::max(metalFraction*cellLength, VSMALL);
        const scalar volumetricHeat = heatFlux/meltThickness;  // [W/m^3]

        // A positive mass flux (evaporation) must *remove* lattice energy.
        // Store the volumetric temperature rate with a sign convention where
        // evaporation yields a negative contribution and condensation positive.
        const scalar rate = volumetricHeat/Cl;
        const scalar relax = 1.0/relaxationTime_;

        phaseChangeRelaxCoeff_[cellI] = relax;
        phaseChangeSource_[cellI] = -rate;
        phaseChangeMassFlux_[cellI] = j_net;
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
    scalar rho1Ref = rho1_.value();
    if (rho1Ref <= SMALL)
    {
        rho1Ref = Foam::gMax(rho1Tmp().internalField());
    }
    forAll(dgdt, cellI)
    {
        const scalar source = phaseChangeSource_[cellI];
        if (Foam::mag(source) <= SMALL)
        {
            dgdt[cellI] = 0.0;
            continue;
        }
        scalar localRate = -(ClVal*source)/(LVal*rho1Ref);
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
        const Foam::scalar minAllowedT = Foam::SMALL;
        if (Foam::min(T) <= minAllowedT)
        {
            WarningInFunction
                << "Received non-positive initial temperature guess."
                << " Clamping to " << minAllowedT << " K to avoid"
                << " thermo inversion failure." << Foam::endl;

            forAll(T, i)
            {
                T[i] = Foam::max(T[i], minAllowedT);
            }
        }
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
    const volScalarField& alpha1Field = alpha1();
    const volScalarField& alpha2Field = alpha2();

    tmp<volScalarField> tCp1 = thermo1_->Cp();
    tmp<volScalarField> tCp2 = thermo2_->Cp();
    tmp<volScalarField> trho1 = thermo1_->rho();
    tmp<volScalarField> trho2 = thermo2_->rho();

    tmp<volScalarField> numerator =
        alpha1Field*trho1()*tCp1()
      + alpha2Field*trho2()*tCp2();

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
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const fvPatchScalarField& alpha1Patch = alpha1().boundaryField()[patchi];
    const fvPatchScalarField& alpha2Patch = alpha2().boundaryField()[patchi];

    tmp<scalarField> tCp1 = thermo1_->Cp(p, T, patchi);
    tmp<scalarField> tCp2 = thermo2_->Cp(p, T, patchi);
    tmp<scalarField> trho1 = thermo1_->rho(patchi);
    tmp<scalarField> trho2 = thermo2_->rho(patchi);

    scalarField result(tCp1().size(), 0.0);

    forAll(result, facei)
    {
        result[facei] = massWeighted
        (
            alpha1Patch[facei],
            trho1()[facei],
            tCp1()[facei],
            alpha2Patch[facei],
            trho2()[facei],
            tCp2()[facei]
        );
    }

    return tmp<scalarField>(new scalarField(std::move(result)));
}
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cp
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

    tmp<scalarField> tCp1 = thermo1_->Cp(p, T, cells);
    tmp<scalarField> tCp2 = thermo2_->Cp(p, T, cells);
    tmp<scalarField> trho1 = thermo1_->rhoEoS(p, T, cells);
    tmp<scalarField> trho2 = thermo2_->rhoEoS(p, T, cells);

    scalarField result(nCells, 0.0);

    forAll(result, i)
    {
        result[i] = massWeighted
        (
            alpha1Field[i],
            trho1()[i],
            tCp1()[i],
            alpha2Field[i],
            trho2()[i],
            tCp2()[i]
        );
    }

    return tmp<scalarField>(new scalarField(std::move(result)));
}
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::Cv() const
{
    const volScalarField& alpha1Field = alpha1();
    const volScalarField& alpha2Field = alpha2();

    tmp<volScalarField> tCv1 = thermo1_->Cv();
    tmp<volScalarField> tCv2 = thermo2_->Cv();
    tmp<volScalarField> trho1 = thermo1_->rho();
    tmp<volScalarField> trho2 = thermo2_->rho();

    tmp<volScalarField> numerator =
        alpha1Field*trho1()*tCv1()
      + alpha2Field*trho2()*tCv2();

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
Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    const fvPatchScalarField& alpha1Patch = alpha1().boundaryField()[patchi];
    const fvPatchScalarField& alpha2Patch = alpha2().boundaryField()[patchi];

    tmp<scalarField> tCv1 = thermo1_->Cv(p, T, patchi);
    tmp<scalarField> tCv2 = thermo2_->Cv(p, T, patchi);
    tmp<scalarField> trho1 = thermo1_->rho(patchi);
    tmp<scalarField> trho2 = thermo2_->rho(patchi);

    scalarField result(tCv1().size(), 0.0);

    forAll(result, facei)
    {
        result[facei] = massWeighted
        (
            alpha1Patch[facei],
            trho1()[facei],
            tCv1()[facei],
            alpha2Patch[facei],
            trho2()[facei],
            tCv2()[facei]
        );
    }

    return tmp<scalarField>(new scalarField(std::move(result)));
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
