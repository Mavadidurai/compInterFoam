/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Copyright (C) 2024
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

Description
    Implementation of the two-temperature model for femtosecond laser-material
    interaction in LIFT processes, providing coupled electron-lattice
    temperature evolution, property management, and diagnostics for LIFT
    simulations.
\*---------------------------------------------------------------------------*/
#include "twoTemperatureModel.H"
#include "fvc.H"
#include "fvm.H"
#include "Pstream.H"
#include "error.H"
#include <cmath>
extern Foam::Switch verbose;
namespace Foam
{
/*---------------------------------------------------------------------------*\
                    twoTemperatureModel Implementation
\*---------------------------------------------------------------------------*/
twoTemperatureModel::twoTemperatureModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField& metalFraction
)
:
    mesh_(mesh),
    dict_(dict),
    metalFraction_(metalFraction),
    ambientTemperature_
    (
        "ambientTemperature",
        dimTemperature,
        dict.lookupOrDefault<scalar>
        (
            "ambientTemperature",
            dict.lookupOrDefault<scalar>("minTe", 300.0)
        )
    ),
    Te_
    (
        IOobject
        (
            "Te",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Te", dimTemperature, 300.0)
    ),
    Tl_
    (
        IOobject
        (
            "Tl",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tl", dimTemperature, 300.0)
    ),
    Ce_
    (
        "Ce",
        dimEnergy/dimVolume/dimTemperature,
        1.5e6
    ),
    Cl_
    (
        "Cl",
        dimEnergy/dimVolume/dimTemperature,
        2.3e6
    ),
    G_
    (
        "G",
        dimEnergy/dimVolume/dimTime/dimTemperature,
        0.0
    ),
    De_
    (
        "De",
        dimLength*dimLength/dimTime,
        1e-4
    ),
    gasMetalExchangeCoeff_
    (
        "gasMetalExchangeCoeff",
        dimEnergy/dimVolume/dimTime/dimTemperature,
        0.0
    ),
    CeFunction_(nullptr),
    GFunction_(nullptr),
    gasMetalExchangeFunction_(nullptr),
    useKapitzaExchange_(false),
    kapitzaZMetal_(0.0),
    kapitzaZGas_(0.0),
    cumulativeLaserEnergy_
    (
        "cumulativeLaserEnergy",
        dimEnergy,
        0.0
    ),
    lastTotalEnergy_
    (
        "lastTotalEnergy",
        dimEnergy,
        0.0
    ),
    lastLoggedEnergy_
    (
        "lastLoggedEnergy",
        dimEnergy,
        0.0
    ),
    energyInitialized_(false),
    loggedEnergyInitialized_(false),
    energyTrackingTimeIndex_(mesh.time().timeIndex()),
    lastElectronSubCycles_(1)
{
    scalar CeLogTe = ambientTemperature_.value();

    if (dict_.found("Ce"))
    {
        if (dict_.isDict("Ce"))
        {
            const dictionary& CeDict = dict_.subDict("Ce");
            CeFunction_.reset(Function1<scalar>::New("Ce", dict_).ptr());

            const scalar minTe = dict_.lookupOrDefault<scalar>("minTe", 300.0);
            const scalar maxTe = dict_.lookupOrDefault<scalar>("maxTe", 4000.0);
            scalar refTe = CeDict.lookupOrDefault<scalar>
            (
                "referenceTemperature",
                CeLogTe
            );

            refTe = Foam::max(minTe, Foam::min(maxTe, refTe));

            const scalar CeAtMin = CeFunction_->value(minTe);
            const scalar CeAtMax = CeFunction_->value(maxTe);
            const scalar CeAtRef = CeFunction_->value(refTe);

            if (CeAtMin <= SMALL || CeAtMax <= SMALL)
            {
                FatalIOErrorInFunction(CeDict)
                    << "Ce Function1 must remain positive in the range ["
                    << minTe << ", " << maxTe << "] K."
                    << exit(FatalIOError);
            }

            if (CeAtRef <= SMALL)
            {
                FatalIOErrorInFunction(CeDict)
                    << "Ce Function1 reference value at Te=" << refTe
                    << " K is non-positive"
                    << exit(FatalIOError);
            }

            Ce_ = dimensionedScalar
            (
                Ce_.name(),
                Ce_.dimensions(),
                CeAtRef
            );

            CeLogTe = refTe;
        }

        else
        {
            dict_.lookup("Ce") >> Ce_;
        }
    }
    else
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'Ce' in two-temperature properties"
            << exit(FatalIOError);
    }

    if (dict_.found("Cl"))
    {
        Cl_ = dict_.lookupOrDefault<dimensionedScalar>("Cl", Cl_);
    }

    if (dict_.found("G"))
    {
        if (dict_.isDict("G"))
        {
            GFunction_.reset(Function1<scalar>::New("G", dict_).ptr());
            const scalar minTemp = dict_.lookupOrDefault<scalar>
            (
                "minTl",
                dict_.lookupOrDefault<scalar>("minTe", ambientTemperature_.value())
            );

            const scalar maxTemp = dict_.lookupOrDefault<scalar>
            (
                "maxTl",
                dict_.lookupOrDefault<scalar>("maxTe", ambientTemperature_.value())
            );

            scalar refTemp = ambientTemperature_.value();
            refTemp = Foam::max(minTemp, Foam::min(maxTemp, refTemp));

            const scalar GAtMin = GFunction_->value(minTemp);
            const scalar GAtMax = GFunction_->value(maxTemp);
            const scalar GAtRef = GFunction_->value(refTemp);

            if (GAtMin <= SMALL || GAtMax <= SMALL || GAtRef <= SMALL)
            {
                FatalIOErrorInFunction(dict_)
                    << "G Function1 must remain positive between "
                    << "minTl/maxTl (" << minTemp << "-" << maxTemp
                    << " K) and at the reference temperature " << refTemp << " K"
                    << exit(FatalIOError);
            }

            G_ = dimensionedScalar(G_.name(), G_.dimensions(), GAtRef);
        }
        else
        {
            dict_.lookup("G") >> G_;
        }
    }
    else
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'G' in two-temperature properties"
            << exit(FatalIOError);
    }

    if (dict_.found("De"))
    {
        De_ = dict_.lookupOrDefault<dimensionedScalar>("De", De_);
    }

    if (dict_.found("gasMetalExchangeCoeff"))
    {
        if (dict_.isDict("gasMetalExchangeCoeff"))
        {
            const dictionary& gasMetalDict =
                dict_.subDict("gasMetalExchangeCoeff");

            word modelType("constant");

            if (gasMetalDict.found("type"))
            {
                gasMetalDict.lookup("type") >> modelType;
            }

            if (modelType == "kapitza" || modelType == "Kapitza")
            {
                useKapitzaExchange_ = true;
                gasMetalExchangeFunction_.clear();
                kapitzaZMetal_ = gasMetalDict.lookupOrDefault<scalar>
                (
                    "Z_metal",
                    2.3e7
                );
                kapitzaZGas_ = gasMetalDict.lookupOrDefault<scalar>
                (
                    "Z_gas",
                    383.0
                );
            }
            else
            {
                gasMetalExchangeFunction_.reset
                (
                    Function1<scalar>::New
                    (
                        "gasMetalExchangeCoeff",
                        dict_
                    ).ptr()
                );
            }
        }
        else
        {
            dict_.lookup("gasMetalExchangeCoeff") >> gasMetalExchangeCoeff_;
        }
    }
    else
    {
        FatalIOErrorInFunction(dict_)
            << "Missing required entry 'gasMetalExchangeCoeff' in"
            << " two-temperature properties"
            << exit(FatalIOError);
    }
    if (!validateParameters())
    {
        FatalErrorInFunction
            << "Invalid model parameters"
            << abort(FatalError);
    }
    // Initialize energy tracking
    updateEnergyTracking();
     const bool master = Pstream::master();
    if (verbose && master)
    {
        Info<< "Two-temperature model initialized:" << nl;

        if (CeFunction_.valid())
        {
            Info<< "  Ce(T) reference @ Te=" << CeLogTe
                << " K = " << Ce_.value() << " J/m³/K" << nl;
        }
        else
        {
            Info<< "  Ce = " << Ce_.value() << " J/m³/K" << nl;
        }

        Info
            << "  Cl = " << Cl_.value() << " J/m³/K" << nl
            << "  G = " << G_.value() << " W/m³/K" << nl
            << "  De = " << De_.value() << " m²/s" << nl
            << "  Gas-metal exchange coeff = "
            << gasMetalExchangeCoeff_.value() << " W/m³/K" << nl
            << "  Ambient temperature = " << ambientTemperature_.value() << " K" << nl
            << "  Metal fraction field = " << metalFraction_.name() << endl;
    }
}
volScalarField& twoTemperatureModel::Te()
{
    if (!validateFields())
    {
        FatalErrorInFunction
            << "Invalid Te/Tl values prior to non-const Te() access"
            << abort(FatalError);
    }
    return Te_;
}
volScalarField& twoTemperatureModel::Tl()
{
    if (!validateFields())
    {
        FatalErrorInFunction
            << "Invalid Te/Tl values prior to non-const Tl() access"
            << abort(FatalError);
    }
    return Tl_;
}
bool twoTemperatureModel::validateParameters() const
{
    // Check dimensions
    if (Te_.dimensions() != dimTemperature ||
        Tl_.dimensions() != dimTemperature)
    {
        FatalErrorInFunction
            << "Invalid temperature dimensions"
            << abort(FatalError);
        return false;
    }
    if (Ce_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        Cl_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        G_.dimensions() != dimEnergy/dimVolume/dimTime/dimTemperature ||
        De_.dimensions() != dimLength*dimLength/dimTime)
    {
        FatalErrorInFunction
            << "Invalid material property dimensions"
            << abort(FatalError);
        return false;
    }
    // Check property values
    if (Ce_.value() <= 0 || Cl_.value() <= 0 || G_.value() <= 0 || De_.value() <= 0)
    {
        FatalErrorInFunction
            << "Non-positive material properties detected"
            << abort(FatalError);
        return false;
    }
	    // Warn if material properties appear outside typical ranges
    if (Ce_.value() < 1e4 || Ce_.value() > 1e8)
    {
        WarningInFunction
            << "Ce (" << Ce_.value() << " J/m³/K) outside typical range [1e4, 1e8]"
            << endl;
    }
    if (Cl_.value() < 1e6 || Cl_.value() > 1e8)
    {
        WarningInFunction
            << "Cl (" << Cl_.value() << " J/m³/K) outside typical range [1e6, 1e8]"
            << endl;
    }
    if (G_.value() < 1e15 || G_.value() > 1e19)
    {
        WarningInFunction
            << "G (" << G_.value() << " W/m³/K) outside typical range [1e15, 1e19]"
            << endl;
    }
    // Validate thermal conductivity dimensions from Ce and De
    if ((Ce_ * De_).dimensions() != dimPower/dimLength/dimTemperature)
    {
        FatalErrorInFunction
            << "Inconsistent thermal conductivity dimensions"
            << abort(FatalError);
        return false;
    }
    // Check temperature values
    forAll(Te_, cellI)
    {
        if (Te_[cellI] < 0 || Tl_[cellI] < 0)
        {
            FatalErrorInFunction
                << "Negative temperature detected at cell " << cellI
                << abort(FatalError);
            return false;
        }
    }
    return true;
}
scalar twoTemperatureModel::metalValidationThreshold() const
{
    const scalar metalFractionFloor =
        dict_.lookupOrDefault<scalar>("metalFractionFloor", 1e-6);
    return dict_.lookupOrDefault<scalar>
    (
        "metalFractionCutoff",
        metalFractionFloor
    );
}
bool twoTemperatureModel::validateFields() const
{
    const scalar metalThreshold = metalValidationThreshold();
    const volScalarField& metalField = metalFraction_;
    auto fieldValid =
        [&](const volScalarField& fld)
    {
        const scalarField& internal = fld.internalField();
        const scalarField& metalInternal = metalField.internalField();
        forAll(internal, cellI)
        {
            if (metalInternal[cellI] < metalThreshold)
            {
                continue;
            }
            const scalar value = internal[cellI];
            if (value <= 0 || !std::isfinite(value))
            {
                return false;
            }
        }
        const volScalarField::Boundary& bf = fld.boundaryField();
        const volScalarField::Boundary& metalBoundary =
            metalField.boundaryField();
        forAll(bf, patchI)
        {
            const fvPatchScalarField& pf = bf[patchI];
            const fvPatchScalarField& metalPatch = metalBoundary[patchI];
            forAll(pf, faceI)
            {
                if (metalPatch[faceI] < metalThreshold)
                {
                    continue;
                }
                const scalar value = pf[faceI];
                if (value <= 0 || !std::isfinite(value))
                {
                    return false;
                }
            }
        }
        return true;
    };
    return fieldValid(Te_) && fieldValid(Tl_);
}

tmp<volScalarField> twoTemperatureModel::metalActiveMask(scalar cutoff) const
{
    const scalar blendWidth =
        dict_.lookupOrDefault<scalar>("metalAmbientBlendWidth", 1e-3);
    const dimensionedScalar cutoffDim("metalCutoff", dimless, cutoff);

    if (blendWidth <= SMALL)
    {
        return pos(metalFraction_ - cutoffDim);
    }

    const scalar lowerBoundValue = Foam::max(cutoff - blendWidth, scalar(0));
    const dimensionedScalar lowerBound
    (
        "metalLowerBound",
        dimless,
        lowerBoundValue
    );

    const dimensionedScalar blendDim
    (
        "metalBlendWidth",
        dimless,
        Foam::max(blendWidth, SMALL)
    );

    tmp<volScalarField> tMask =
        Foam::min
        (
            dimensionedScalar("one", dimless, 1.0),
            Foam::max
            (
                dimensionedScalar("zero", dimless, 0.0),
                (metalFraction_ - lowerBound)/blendDim
            )
        );

    return tMask;
}
tmp<volScalarField> twoTemperatureModel::clampedMetalFraction() const
{
    const dimensionedScalar metalZero("metalFractionZero", dimless, 0.0);
    const dimensionedScalar metalCeiling("metalFractionCeiling", dimless, 1.0);

    return Foam::min(Foam::max(metalFraction_, metalZero), metalCeiling);
}

void twoTemperatureModel::applyTemperatureBounds
(
    const volScalarField& activeMask,
    const dimensionedScalar& minTe,
    const dimensionedScalar& maxTe,
    const dimensionedScalar& minTl,
    const dimensionedScalar& maxTl,
    const dimensionedScalar& ambient
)
{
    tmp<volScalarField> tBoundTe = Foam::max(Foam::min(Te_, maxTe), minTe);
    tmp<volScalarField> tBoundTl = Foam::max(Foam::min(Tl_, maxTl), minTl);
    const volScalarField& boundTe = tBoundTe();
    const volScalarField& boundTl = tBoundTl();

    const dimensionedScalar activeThreshold
    (
        "activeThreshold",
        dimless,
        VSMALL
    );

    tmp<volScalarField> tBinaryActive = pos(activeMask - activeThreshold);
    const volScalarField& binaryActive = tBinaryActive();

    tmp<volScalarField> tInactiveMask = scalar(1) - binaryActive;
    const volScalarField& inactiveMask = tInactiveMask();

    Te_ = binaryActive*boundTe + inactiveMask*ambient;
    Tl_ = binaryActive*boundTl + inactiveMask*ambient;

    Te_.correctBoundaryConditions();
    Tl_.correctBoundaryConditions();
}

void twoTemperatureModel::solveLatticeEquation
(
    const volScalarField& metalEff,
    const volScalarField& metalPhysical,
    const volScalarField& G,
    const volScalarField& klField,
    const volScalarField& phaseChangeRelaxCoeff,
    const volScalarField& phaseChangeSource,
    const volScalarField& TlOld,
    const volScalarField& gasMetalHeatFlux,
    const dimensionedScalar& dtSub,
    const volScalarField& TlPrev
)
{
    tmp<volScalarField> tCap = metalEff*Cl_;
    const volScalarField& capacity = tCap();

    fvScalarMatrix TlEqn
    (
        fvm::Sp(capacity/dtSub, Tl_)
      - fvm::laplacian(metalEff*klField, Tl_)
      + fvm::Sp(metalEff*G, Tl_)
      + fvm::Sp(metalEff*Cl_*phaseChangeRelaxCoeff, Tl_)
    ==
        (capacity/dtSub)*TlPrev
      + metalPhysical*G*Te_
      + metalPhysical*Cl_*phaseChangeSource
      + metalPhysical*Cl_*phaseChangeRelaxCoeff*TlOld
      - metalPhysical*gasMetalHeatFlux
    );

    TlEqn.relax();
    TlEqn.solve(mesh_.solver("Tl"));
    Tl_.correctBoundaryConditions();
}

void twoTemperatureModel::solveElectronEquation
(
    const volScalarField& metalEff,
    const volScalarField& metalPhysical,
    const volScalarField& Ce,
    const volScalarField& ke,
    const volScalarField& G,
    const volScalarField& laserSource,
    const dimensionedScalar& dtSub,
    const volScalarField& TePrev
)
{
    tmp<volScalarField> tCap = metalEff*Ce;
    const volScalarField& capacity = tCap();

    fvScalarMatrix TeEqn
    (
        fvm::Sp(capacity/dtSub, Te_)
      - fvm::laplacian(metalEff*ke, Te_)
      + fvm::Sp(metalEff*G, Te_)
    ==
        (capacity/dtSub)*TePrev
      + metalPhysical*laserSource
      + metalPhysical*G*Tl_
    );

    TeEqn.relax();
    TeEqn.solve(mesh_.solver("Te"));
    Te_.correctBoundaryConditions();
}

dimensionedScalar twoTemperatureModel::couplingEnergy
(
    const volScalarField& gasMetalHeatFlux,
    const dimensionedScalar& dtDim
) const
{
    tmp<volScalarField> tMetal = clampedMetalFraction();
    const volScalarField& metalEff = tMetal();
    return fvc::domainIntegrate(metalEff*gasMetalHeatFlux)*dtDim;
    
}

dimensionedScalar twoTemperatureModel::currentTotalEnergy() const
{
    tmp<volScalarField> tCe = electronHeatCapacity();
    const volScalarField& CeField = tCe();
    tmp<volScalarField> tMetal = clampedMetalFraction();
    const volScalarField& metalEff = tMetal();
    return fvc::domainIntegrate(metalEff*CeField*Te_ + metalEff*Cl_*Tl_);
}

void twoTemperatureModel::writeEnergyDiagnostics
(
    const dimensionedScalar& laserEnergy,
    const dimensionedScalar& phaseChangeEnergy,
    const dimensionedScalar& couplingEnergy,
    const dimensionedScalar& electronEnergyBefore,
    const dimensionedScalar& latticeEnergyBefore,
    const Switch& energyDiagnostics
) const
{
    if (!energyDiagnostics)
    {
        return;
    }

    tmp<volScalarField> tCeAfter = electronHeatCapacity();
    const volScalarField& CeAfter = tCeAfter();
    tmp<volScalarField> tMetal = clampedMetalFraction();
    const volScalarField& metalEff = tMetal();
    const dimensionedScalar electronEnergyAfter =
        fvc::domainIntegrate(metalEff*CeAfter*Te_);
    const dimensionedScalar latticeEnergyAfter =
        fvc::domainIntegrate(metalEff*Cl_*Tl_);

    cumulativeLaserEnergy_ += laserEnergy;

    Info<< "Energy diagnostics:" << nl
        << "  Electron energy before: "
        << electronEnergyBefore.value() << " J" << nl
        << "  Lattice energy before: "
        << latticeEnergyBefore.value() << " J" << nl
        << "  Laser energy input: " << laserEnergy.value() << " J" << nl
        << "  Phase-change energy input: "
        << phaseChangeEnergy.value() << " J" << nl
        << "  Gas-metal coupling energy loss: "
        << couplingEnergy.value() << " J" << nl
        << "  Electron energy after: "
        << electronEnergyAfter.value() << " J" << nl
        << "  Lattice energy after: "
        << latticeEnergyAfter.value() << " J" << nl
        << "  Total metal energy: "
        << (electronEnergyAfter + latticeEnergyAfter).value() << " J" << nl
        << "  Cumulative laser energy: "
        << cumulativeLaserEnergy_.value() << " J" << endl;
        
    const scalar laserEnergyValue = laserEnergy.value();
    const scalar phaseChangeEnergyValue = phaseChangeEnergy.value();

    Info<< "  Phase-change/Laser energy ratio: ";
    if (mag(laserEnergyValue) > SMALL)
    {
        const scalar ratio = phaseChangeEnergyValue/laserEnergyValue;
        Info<< ratio << nl;

        if (ratio > 2.0)
        {
            WarningInFunction
                << "Phase-change energy input exceeds twice the laser energy input." << nl
                << "  Phase-change energy: " << phaseChangeEnergyValue << " J" << nl
                << "  Laser energy: " << laserEnergyValue << " J" << nl
                << "Verify phase-change settings and source implementation." << endl;
        }
    }
    else if (mag(phaseChangeEnergyValue) > SMALL)
    {
        Info<< "undefined (laser energy ≈ 0)" << nl;
        WarningInFunction
            << "Non-negligible phase-change energy reported while laser energy is near zero." << nl
            << "  Phase-change energy: " << phaseChangeEnergyValue << " J" << endl;
        return;
    }
    else
    {
        Info<< "undefined (laser energy ≈ 0)" << nl;
        return;
    }

    Info<< flush;
}

void twoTemperatureModel::writeSolveStatistics
(
    const scalar residual,
    const Switch& statsSwitch
) const
{
    if (!statsSwitch)
    {
        return;
    }

    const dimensionedScalar minTe = min(Te_);
    const dimensionedScalar maxTe = max(Te_);
    const dimensionedScalar minTl = min(Tl_);
    const dimensionedScalar maxTl = max(Tl_);

    tmp<volScalarField> tTempDiff = mag(Te_ - Tl_);
    const volScalarField& tempDiff = tTempDiff();

    Info<< "Two-temperature solve:" << nl
        << "  Te range: " << minTe.value() << " - " << maxTe.value() << " K" << nl
        << "  Tl range: " << minTl.value() << " - " << maxTl.value() << " K" << nl
        << "  Max temperature difference: "
        << max(tempDiff).value() << " K" << nl
        << "  Mean Te: " << gAverage(Te_) << " K" << nl
        << "  Mean Tl: " << gAverage(Tl_) << " K" << nl
        << "  Coupling residual: " << residual << endl;
}

bool twoTemperatureModel::checkEnergyConservation
(
    const dimensionedScalar& previousEnergy,
    const dimensionedScalar& expectedEnergyChange,
    const dimensionedScalar& currentEnergy,
    scalar& relativeError
) const
{
    if (!energyInitialized_)
    {
        relativeError = 0.0;
        return true;
    }

    const dimensionedScalar expectedEnergy =
        previousEnergy + expectedEnergyChange;
    relativeError = mag
    (
        (currentEnergy.value() - expectedEnergy.value())/
        (mag(expectedEnergy.value()) + SMALL)
    );

    const scalar energyTol =
        dict_.lookupOrDefault<scalar>
        (
            "energyTol",
            dict_.lookupOrDefault<scalar>("energyTolerance", 0.01)
        );
    return relativeError < energyTol;
}

void twoTemperatureModel::reportEnergyViolation
(
    const dimensionedScalar& previousEnergy,
    const dimensionedScalar& totalEnergyInput,
    const dimensionedScalar& currentEnergy,
    scalar relativeError,
    scalar dt,
    const Switch& auditSwitch
) const
{
    const dimensionedScalar expectedEnergy = previousEnergy + totalEnergyInput;

    if (!auditSwitch)
    {
        WarningInFunction
            << "Energy conservation violation detected. Relative error = "
            << relativeError*100 << " %" << endl;
        return;
    }

    const scalar tauE = Ce_.value()/G_.value();
    const scalar tauL = Cl_.value()/G_.value();

    WarningInFunction
        << "Energy conservation violation detected" << nl
        << "Error = " << relativeError * 100 << " %" << nl
        << "Previous energy: " << previousEnergy.value() << " J" << nl
        << "Net source energy: " << totalEnergyInput.value() << " J" << nl
        << "Expected energy: " << expectedEnergy.value() << " J" << nl
        << "Current energy: " << currentEnergy.value() << " J" << nl
        << "Ce = " << Ce_.value() << " J/m^3/K" << nl
        << "Cl = " << Cl_.value() << " J/m^3/K" << nl
        << "G  = " << G_.value()  << " W/m^3/K" << nl
        << "deltaT = " << dt << " s" << nl
        << "Characteristic times: Ce/G = " << tauE
        << " s, Cl/G = " << tauL << " s" << nl
        << "Suggested max deltaT: " << 0.2*Foam::min(tauE, tauL)
        << " s" << endl;
}

void twoTemperatureModel::updateEnergyTracking() const
{
    updateEnergyTracking(currentTotalEnergy());
}

void twoTemperatureModel::updateEnergyTracking
(
    const dimensionedScalar& currentEnergy
) const
{

    const label currentTimeIndex = mesh_.time().timeIndex();

    if (!energyInitialized_ || currentTimeIndex < energyTrackingTimeIndex_)
    {
        cumulativeLaserEnergy_ = dimensionedScalar
        (
            cumulativeLaserEnergy_.name(),
            cumulativeLaserEnergy_.dimensions(),
            0.0
        );
    }
    lastTotalEnergy_ = currentEnergy;
    energyInitialized_ = true;
    loggedEnergyInitialized_ = false;
    energyTrackingTimeIndex_ = currentTimeIndex;   
}

label twoTemperatureModel::electronSubCycleCount
(
    const dimensionedScalar& dt
) const
{
    label planned = dict_.lookupOrDefault<label>("electronSubCycles", 1);

    const scalar dtValue = Foam::max(dt.value(), VSMALL);

    if (dict_.found("maxElectronDeltaT"))
    {
        const dimensionedScalar defaultMaxElectronDt
        (
            "maxElectronDeltaT",
            dimTime,
            dtValue
        );

        dimensionedScalar maxElectronDtDim =
            dict_.lookupOrDefault<dimensionedScalar>
            (
                "maxElectronDeltaT",
                defaultMaxElectronDt
            );

        if (maxElectronDtDim.dimensions() == dimTime)
        {
            // Already dimensionally consistent
        }
        else if (maxElectronDtDim.dimensions() == dimless)
        {
            maxElectronDtDim = dimensionedScalar
            (
                maxElectronDtDim.name(),
                dimTime,
                maxElectronDtDim.value()
            );
        }
        else
        {
            FatalIOErrorInFunction(dict_)
                << "Entry 'maxElectronDeltaT' must be specified either "
                << "with time dimensions or as a plain scalar." << nl
                << "    dimensions provided: "
                << maxElectronDtDim.dimensions() << nl
                << "    expected: " << dimTime
                << " or " << dimless
                << exit(FatalIOError);
        }


        const scalar maxElectronDt = maxElectronDtDim.value();

        if (maxElectronDt > SMALL)
        {
            const scalar ratio = dtValue/maxElectronDt;
            planned = Foam::max
            (
                planned,
                label(std::ceil(ratio))
            );
        }
    }

    const label minCycles =
        dict_.lookupOrDefault<label>("minElectronSubCycles", 1);
    planned = Foam::max(planned, Foam::max(minCycles, label(1)));

    const label maxCycles =
        dict_.lookupOrDefault<label>("maxElectronSubCycles", 0);

    if (maxCycles > 0)
    {
        planned = Foam::min(planned, maxCycles);
    }

    return Foam::max(planned, label(1));
}

label twoTemperatureModel::activeElectronSubCycles
(
    const dimensionedScalar& dt
)
{
    const label planned = electronSubCycleCount(dt);
    lastElectronSubCycles_ = planned;
    return planned;
}

label twoTemperatureModel::plannedElectronSubCycles
(
    const dimensionedScalar& dt
) const
{
    return electronSubCycleCount(dt);
}
void twoTemperatureModel::solve
(
    const volScalarField& laserSource,
    const volScalarField& phaseChangeSource,
    const volScalarField& phaseChangeRelaxCoeff,
    const volScalarField& gasMetalHeatFlux
)
{
    if (!validateFields())
    {
        FatalErrorInFunction
            << "Invalid field values before solve"
            << abort(FatalError);
    }

    const volScalarField& metal = metalFraction_;

    const scalar ambient = ambientTemperature_.value();
    const dimensionedScalar ambientDim
    (
        "ambientTemperature",
        dimTemperature,
        ambient
    );
    const scalar metalFractionFloor =
        dict_.lookupOrDefault<scalar>("metalFractionFloor", 1e-6);
    const scalar metalCutoff =
        dict_.lookupOrDefault<scalar>
        (
            "metalFractionCutoff",
            metalFractionFloor
        );
    const dimensionedScalar metalFloor
    (
        "metalFractionFloor",
        dimless,
        Foam::max(metalFractionFloor, VSMALL)
    );
    const dimensionedScalar metalCeiling("metalFractionCeiling", dimless, 1.0);
    tmp<volScalarField> tMetalEff =
        Foam::min(Foam::max(metal, metalFloor), metalCeiling);
    const volScalarField& metalEff = tMetalEff();
    const dimensionedScalar metalZero("metalFractionZero", dimless, 0.0);
    tmp<volScalarField> tMetalPhysical =
        Foam::min(Foam::max(metal, metalZero), metalCeiling);
    const volScalarField& metalPhysical = tMetalPhysical();

    const dimensionedScalar minTe
    (
        "minTe",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("minTe", 300.0)
    );
    const dimensionedScalar maxTe
    (
        "maxTe",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("maxTe", 3500.0)
    );
    const dimensionedScalar minTl
    (
        "minTl",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("minTl", minTe.value())
    );
    const dimensionedScalar maxTl
    (
        "maxTl",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("maxTl", maxTe.value())
    );

    const Switch energyDiagnostics
    (
        dict_.lookupOrDefault<Switch>("energyDiagnostics", false)
    );
    const Switch temperatureDiagnostics
    (
        dict_.lookupOrDefault<Switch>("temperatureDiagnostics", verbose)
    );
    const Switch energyAudit
    (
        dict_.lookupOrDefault<Switch>("energyAudit", verbose)
    );

    dimensionedScalar previousEnergy("previousEnergy", dimEnergy, 0.0);
    if (energyInitialized_)
    {
        previousEnergy = lastTotalEnergy_;
    }
    else
    {
        previousEnergy = currentTotalEnergy();
    }
    tmp<volScalarField> tActiveMask = metalActiveMask(metalCutoff);
    const volScalarField& activeMask = tActiveMask();
    applyTemperatureBounds(activeMask, minTe, maxTe, minTl, maxTl, ambientDim);

    // Recompute the baseline energy after enforcing bounds so that
    // redistribution caused by updated metal fractions does not appear
    // as an artificial conservation error.
    previousEnergy = currentTotalEnergy();

    const dimensionedScalar dtDim = mesh_.time().deltaT();
    // Cache the micro-step count for diagnostics; callers can query the
    // const plannedElectronSubCycles() helper without altering this record.
    const label nElectronSubCycles =
        Foam::max(activeElectronSubCycles(dtDim), label(1));
    const dimensionedScalar dtSub = dtDim/scalar(nElectronSubCycles);

    tmp<volScalarField> tTePrev
    (
        new volScalarField
        (
            IOobject
            (
                "TePrevSub",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("TePrevSub", dimTemperature, 0.0)
        )
    );
    volScalarField& TePrev = tTePrev.ref();

    tmp<volScalarField> tTlPrev
    (
        new volScalarField
        (
            IOobject
            (
                "TlPrevSub",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("TlPrevSub", dimTemperature, 0.0)
        )
    );
    volScalarField& TlPrev = tTlPrev.ref();
    const dimensionedScalar laserEnergy =
        fvc::domainIntegrate(metalPhysical*laserSource)*dtDim;
    const dimensionedScalar phaseChangeEnergy =
        fvc::domainIntegrate(metalPhysical*Cl_*phaseChangeSource)*dtDim;
    const dimensionedScalar couplingEnergyValue =
        couplingEnergy(gasMetalHeatFlux, dtDim);
    dimensionedScalar clampEnergyCorrection
    (
        "clampEnergyCorrection",
        dimEnergy,
        0.0
    );
    const dimensionedScalar totalEnergyInput =
        laserEnergy + phaseChangeEnergy - couplingEnergyValue;

    dimensionedScalar electronEnergyBefore("electronEnergyBefore", dimEnergy, 0.0);
    dimensionedScalar latticeEnergyBefore("latticeEnergyBefore", dimEnergy, 0.0);

    if (energyDiagnostics)
    {
        tmp<volScalarField> tCeInitial = electronHeatCapacity();
        const volScalarField& CeInitial = tCeInitial();
        electronEnergyBefore =
            fvc::domainIntegrate(metalPhysical*CeInitial*Te_);
        latticeEnergyBefore =
            fvc::domainIntegrate(metalPhysical*Cl_*Tl_);
    }

    const label nInnerSweeps =
        dict_.lookupOrDefault<label>("nInnerCouplingSweeps", 1);
    const scalar residualTol =
        dict_.lookupOrDefault<scalar>
        (
            "innerCouplingResidualTol",
            dict_.lookupOrDefault<scalar>
            (
                "innerCouplingResidualTolerance",
                dict_.lookupOrDefault<scalar>
                (
                    "innerCouplingReductionTol",
                    dict_.lookupOrDefault<scalar>
                    (
                        "innerCouplingReductionTolerance",
                        -GREAT
                    )
                )
            )
        );

    scalar residual = gMax(mag(Te_ - Tl_)().internalField());
    const volScalarField& TlOld = Tl_.oldTime();

    for (label sweep = 0; sweep < nInnerSweeps; ++sweep)
    {
        TePrev = Te_.oldTime();
        TlPrev = TlOld;

        for (label sub = 0; sub < nElectronSubCycles; ++sub)
        {
            tmp<volScalarField> tG = electronPhononCoupling();
            const volScalarField& G = tG();

            tmp<volScalarField> tkl = kl();
            const volScalarField& klField = tkl();

            solveLatticeEquation
            (
                metalEff,
                metalPhysical,
                G,
                klField,
                phaseChangeRelaxCoeff,
                phaseChangeSource,
                TlOld,
                gasMetalHeatFlux,
                dtSub,
                TlPrev
            );
            tmp<volScalarField> tTlClamped =
                Foam::max(Foam::min(Tl_, maxTl), minTl);
            const volScalarField& TlClamped = tTlClamped();
            tmp<volScalarField> tTlDelta = Tl_ - TlClamped;
            const volScalarField& TlDelta = tTlDelta();
            clampEnergyCorrection +=
                fvc::domainIntegrate(metalPhysical*Cl_*TlDelta);
            Tl_ = TlClamped;
            Tl_.correctBoundaryConditions();
            TlPrev = Tl_;

            tmp<volScalarField> tke = electronThermalConductivity();
            const volScalarField& keField = tke();

            tmp<volScalarField> tCe = electronHeatCapacity();
            const volScalarField& CeField = tCe();

            solveElectronEquation
            (
                metalEff,
                metalPhysical,
                CeField,
                keField,
                G,
                laserSource,
                dtSub,
                TePrev
            );
            tmp<volScalarField> tCeClamp = electronHeatCapacity();
            const volScalarField& CeClamp = tCeClamp();
            tmp<volScalarField> tTeClamped =
                Foam::max(Foam::min(Te_, maxTe), minTe);
            const volScalarField& TeClamped = tTeClamped();
            tmp<volScalarField> tTeDelta = Te_ - TeClamped;
            const volScalarField& TeDelta = tTeDelta();
            clampEnergyCorrection +=
                fvc::domainIntegrate(metalPhysical*(CeClamp*TeDelta));
            Te_ = TeClamped;
            Te_.correctBoundaryConditions();
            TePrev = Te_;
        }

        residual = gMax(mag(Te_ - Tl_)().internalField());

        if (residualTol > 0 && residual < residualTol)
        {
            if (temperatureDiagnostics)
            {
                Info<< "Inner two-temperature coupling converged after "
                    << sweep + 1 << " sweeps with residual "
                    << residual << endl;
            }
            break;
        }
    }
 // Energy balance diagnostics
    if (verbose && Pstream::master())
    {
        tmp<volScalarField> tCeFinal = electronHeatCapacity();
        //const volScalarField& CeFinal = tCeFinal();

        tmp<volScalarField> tGFinal = electronPhononCoupling();
        const volScalarField& GFinal = tGFinal();

        const dimensionedScalar laserPower =
            fvc::domainIntegrate(metalPhysical*laserSource);
        const dimensionedScalar eCoupling =
            fvc::domainIntegrate(metalPhysical*GFinal*Te_);
        const dimensionedScalar lCoupling =
            fvc::domainIntegrate(metalPhysical*GFinal*Tl_);
        const dimensionedScalar gasLoss =
            fvc::domainIntegrate(metalPhysical*gasMetalHeatFlux);
        const dimensionedScalar metalVolume =
            fvc::domainIntegrate(metalPhysical);
        const scalar dtValue = dtDim.value();
        scalar clampPower = 0.0;
        if (dtValue > SMALL)
        {
            clampPower = clampEnergyCorrection.value()/dtValue;
        }

        Info<< "══════ ENERGY BALANCE ══════" << nl
            << "Metal volume: " << metalVolume.value()*1e18 << " µm³" << nl
            << "Power terms [W]:" << nl
            << "  Laser input:        " << laserPower.value() << nl
            << "  e→l coupling:       " << (eCoupling - lCoupling).value() << nl
            << "  Gas coupling loss:  " << gasLoss.value() << nl
            << "  Clamp correction:   " << clampPower << nl
            << "Net power [W]:" << nl
            << "  Into electrons:     " << (laserPower - eCoupling + lCoupling).value() << nl
            << "  Into lattice:       " << (eCoupling - lCoupling - gasLoss).value() << nl
            << "Temperatures:" << nl
            << "  max(Te): " << gMax(Te_) << " K" << nl
            << "  max(Tl): " << gMax(Tl_) << " K" << nl
            << "════════════════════════════" << endl;
    }

    const dimensionedScalar energyBeforeFinalClamp = currentTotalEnergy();
    applyTemperatureBounds(activeMask, minTe, maxTe, minTl, maxTl, ambientDim);
    const dimensionedScalar energyAfterFinalClamp = currentTotalEnergy();
    clampEnergyCorrection += energyBeforeFinalClamp - energyAfterFinalClamp;
    residual = gMax(mag(Te_ - Tl_)().internalField());

    tmp<volScalarField> tCeFinal = electronHeatCapacity();
    const volScalarField& CeFinal = tCeFinal();
    const dimensionedScalar currentEnergy =
        fvc::domainIntegrate(metalPhysical*(CeFinal*Te_ + Cl_*Tl_));

    scalar energyError = 0.0;
    const dimensionedScalar netEnergyInput =
        totalEnergyInput - clampEnergyCorrection;
    if
    (
        !checkEnergyConservation
        (
            previousEnergy,
            netEnergyInput,
            currentEnergy,
            energyError
        )
    )
    {
        reportEnergyViolation
        (
            previousEnergy,
            netEnergyInput,
            currentEnergy,
            energyError,
            dtDim.value(),
            energyAudit
        );
    }

    const dimensionedScalar previousLoggedEnergy =
        energyInitialized_ ? lastTotalEnergy_ : currentEnergy;

    updateEnergyTracking(currentEnergy);

    lastLoggedEnergy_ = previousLoggedEnergy;
    loggedEnergyInitialized_ = true;

    if (energyDiagnostics && mag(clampEnergyCorrection.value()) > SMALL)
    {
        Info<< "  Temperature clamp energy correction: "
            << clampEnergyCorrection.value() << " J" << endl;
    }
    writeEnergyDiagnostics
    (
        laserEnergy,
        phaseChangeEnergy,
        couplingEnergyValue,
        electronEnergyBefore,
        latticeEnergyBefore,
        energyDiagnostics
    );

    writeSolveStatistics(residual, temperatureDiagnostics);
}

tmp<volScalarField> twoTemperatureModel::electronThermalConductivity() const
{
    tmp<volScalarField> ke
    (
        new volScalarField
        (
            IOobject
            (
                "ke",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("ke", dimPower/dimLength/dimTemperature, 1.0)
        )
    );
    // Get primitive field references for efficient looping
    scalarField& kei = ke.ref().primitiveFieldRef();
    tmp<volScalarField> tCe = electronHeatCapacity();
    const volScalarField& CeField = tCe();
    const scalarField& Cei = CeField.primitiveField();
    forAll(kei, cellI)
    {
        kei[cellI] = (Cei[cellI] * De_).value();
    }
    ke.ref().correctBoundaryConditions();
    return ke;
}
tmp<volScalarField> twoTemperatureModel::electronHeatCapacity() const
{
    // Create temporary field for electron heat capacity
    tmp<volScalarField> tCe
    (
        new volScalarField
        (
            IOobject
            (
                "Ce",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                 IOobject::NO_WRITE,
                false
            ),
            mesh_,
            Ce_
        )
    );
    volScalarField& CeField = tCe.ref();
    if (CeFunction_.valid())
    {
        forAll(CeField, cellI)
        {
            const scalar Te = Te_[cellI];
            const scalar gamma = 630.0;
            const scalar CeValue = Foam::max(gamma*Te, scalar(1e4));
            CeField[cellI] = CeValue;
        }
    }
    else
    {
        const scalar CeConst = Foam::max(Ce_.value(), scalar(1e4));
        CeField = dimensionedScalar(Ce_.name(), Ce_.dimensions(), CeConst);
    }

    CeField.correctBoundaryConditions();

    return tCe;
}
tmp<volScalarField> twoTemperatureModel::electronPhononCoupling() const
{
    // Create temporary field for electron-phonon coupling
    tmp<volScalarField> tG
    (
        new volScalarField
        (
            IOobject
            (
                "G",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            G_
        )
    ); // Electron-phonon coupling field
    volScalarField& GField = tG.ref(); // Reference for modification
    if (GFunction_.valid())
    {
        forAll(GField, cellI)
        {
            const scalar Te = Te_[cellI];
            const scalar Tl = Tl_[cellI];
            const scalar TlSafe = Foam::max(Tl, scalar(1.0));

            const scalar baseG = GFunction_->value(TlSafe);

            const scalar TRatio = Foam::max(Te/TlSafe, scalar(1.0));
            scalar GValue = baseG * Foam::sqrt(TRatio);

            GField[cellI] = Foam::min(GValue, scalar(1e19));
        }
    }
    else
    {
        const scalar baseG = G_.value();
        forAll(GField, cellI)
        {
            const scalar Te = Te_[cellI];
            const scalar Tl = Tl_[cellI];
            const scalar TlSafe = Foam::max(Tl, scalar(1.0));

            const scalar TRatio = Foam::max(Te/TlSafe, scalar(1.0));
            scalar GValue = baseG * Foam::sqrt(TRatio);

            GField[cellI] = Foam::min(GValue, scalar(1e19));
        }
    }

    GField.correctBoundaryConditions();

    return tG;
}
tmp<volScalarField> twoTemperatureModel::gasMetalExchangeCoeffField() const
{
    tmp<volScalarField> tCoeff
    (
        new volScalarField
        (
            IOobject
            (
                "gasMetalExchangeCoeffField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            gasMetalExchangeCoeff_
        )
    );
    
    volScalarField& coeff = tCoeff.ref();

    if (useKapitzaExchange_)
    {
        // Gas-metal Kapitza conductance based on acoustic mismatch
        const volScalarField& Tl = Tl_;

        const scalar zTi = kapitzaZMetal_;
        const scalar zAr = kapitzaZGas_;
        const scalar zRatio = zAr/(zTi + VSMALL);
        const scalar tau = 4.0*zRatio/sqr(1.0 + zRatio);

        const scalar kB = 1.38e-23;
        const scalar hBar = 1.055e-34;
        const scalar pi = 3.141592653589793;
        const scalar prefactor = (kB*kB)/(6.0*sqr(pi)*Foam::pow3(hBar));
        const scalar delta = 1e-9;  // m

        forAll(coeff, cellI)
        {
            const scalar T = Foam::max(Tl[cellI], scalar(0));

            const scalar hK = prefactor*tau*T*T*T;
            const scalar hVol = hK/delta;

            // The interfacial conductance must remain active even when the
            // volume fraction has been clipped to 0 or 1 by MULES.  The
            // additional 4α(1-α) weighting that was introduced earlier drove
            // the coefficient to zero in precisely those sharply resolved
            // interface cells, which is why the runtime log reported a
            // vanishing "mean exchange coeff" and "max |q_gm|".  With the
            // conductance quenched, the gas never absorbed energy from the
            // hot lattice and the simulation ran away thermally.  Removing the
            // spurious weight restores the expected Kapitza exchange whenever
            // both phases occupy the cell, leaving the masking in TEqn.H to
            // gate the coupling elsewhere.

            scalar val = hVol;
            val = Foam::min(val, scalar(1e14));
            val = Foam::max(val, scalar(0));

            coeff[cellI] = val;
        }
    }

    coeff.correctBoundaryConditions();
    return tCoeff;
}
bool twoTemperatureModel::valid() const
{
    if (!validateParameters())
    {
        return false;
    }

    if (!validateFields())
    {
        return false;
    }
    // Check field dimensions
    if (Te_.dimensions() != dimTemperature ||
        Tl_.dimensions() != dimTemperature ||
        Ce_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        Cl_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        G_.dimensions() != dimEnergy/dimVolume/dimTime/dimTemperature ||
        De_.dimensions() != dimLength*dimLength/dimTime)
    {
        return false;
    }
    return true;
}
void twoTemperatureModel::write() const
{
    const bool master = Pstream::master();
    if (verbose && master)
    {
        Info<< "Two-temperature model:" << nl
            << "Parameters:" << nl
            << "  Ce = " << Ce_.value() << " J/m³/K" << nl
            << "  Cl = " << Cl_.value() << " J/m³/K" << nl
            << "  G = " << G_.value() << " W/m³/K" << nl
            << "Field statistics:" << nl
            << "  Te range: " << min(Te_).value() << " - " << max(Te_).value() << " K" << nl
            << "  Tl range: " << min(Tl_).value() << " - " << max(Tl_).value() << " K" << nl
            << "  Mean Te: " << gAverage(Te_) << " K" << nl
            << "  Mean Tl: " << gAverage(Tl_) << " K" << nl;
    }
    if (energyInitialized_ && loggedEnergyInitialized_)
    {
        tmp<volScalarField> tCeW = electronHeatCapacity();
        const volScalarField& CeW = tCeW();
        tmp<volScalarField> tMetal = clampedMetalFraction();
        const volScalarField& metalEff = tMetal();
        dimensionedScalar currentEnergy =
            fvc::domainIntegrate(metalEff*(CeW*Te_ + Cl_*Tl_));
        scalar energyError = mag
        (
            (currentEnergy.value() - lastLoggedEnergy_.value())/
            (mag(lastLoggedEnergy_.value()) + SMALL)
        );
        if (verbose && master)
        {
            Info<< "Energy conservation:" << nl
                << "  Current total energy: " << currentEnergy.value() << " J" << nl
                << "  Energy error: " << energyError * 100 << " %" << endl;
        }
        lastLoggedEnergy_ = currentEnergy;
    }
}
tmp<volScalarField> Foam::twoTemperatureModel::ke() const
{
    return electronThermalConductivity(); // Electronic thermal conductivity from constant diffusivity
}
tmp<volScalarField> Foam::twoTemperatureModel::kl() const
{
    const scalar klHighTThreshold =
        dict_.lookupOrDefault<scalar>("klHighTThreshold", 1000.0);
    const scalar klExponent =
        dict_.lookupOrDefault<scalar>("klExponent", 0.5);    
    tmp<volScalarField> tkl
    (
        new volScalarField
        (
            IOobject
            (
                "kl",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "kl",
                dimPower/dimLength/dimTemperature,
                Cl_.value()*De_.value()
            )
        )
    );
    volScalarField& kl = tkl.ref();
    // Calculate temperature-dependent lattice thermal conductivity
    // Using Drude model with phonon contribution
    forAll(mesh_.C(), cellI)
    {
        scalar Tl = Tl_[cellI];
        // Apply temperature dependent correction
        if (Tl > klHighTThreshold && klExponent != 0.0)
        {
            const scalar ratio =
                Foam::max(klHighTThreshold, VSMALL)/Foam::max(Tl, VSMALL);
            kl[cellI] *= pow(ratio, klExponent);  // High temperature correction
        }
    }
    return tkl;
}
} // End namespace Foam
