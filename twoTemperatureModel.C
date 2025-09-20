/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Implementation of the two-temperature model for femtosecond laser-material
    interaction in LIFT process. Handles:
    - Temperature field initialization and evolution
    - Material property calculations
    - Energy conservation tracking
    - Coupled electron-lattice temperature solution
    The model uses temperature-dependent material properties and includes
    electron-phonon coupling for accurate simulation of ultrafast laser heating.
\*---------------------------------------------------------------------------*/
#include "twoTemperatureModel.H"
#include "fvc.H"
#include "fvm.H"
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
        210.0
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
        5e17
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
        5e17
    ),
    CeFunction_(nullptr),
    GFunction_(nullptr),
    gasMetalExchangeFunction_(nullptr),
    lastTotalEnergy_
    (
        "lastTotalEnergy",
        dimEnergy,
        0.0
    ),
    energyInitialized_(false)
{
if (dict.found("Ce"))
    {
        if (dict.isDict("Ce"))
        {
            CeFunction_.reset(Function1<scalar>::New("Ce", dict.subDict("Ce")).ptr());
        }
        else
        {
            Ce_ = dict.lookupOrDefault<dimensionedScalar>("Ce", Ce_);
        }
    }
    if (dict.found("Cl"))
    {
        Cl_ = dict.lookupOrDefault<dimensionedScalar>("Cl", Cl_);
    }
    if (dict.found("G"))
    {
        if (dict.isDict("G"))
        {
            GFunction_.reset(Function1<scalar>::New("G", dict.subDict("G")).ptr());
        }
        else
        {
            G_ = dict.lookupOrDefault<dimensionedScalar>("G", G_);
        }
    }
    if (dict.found("De"))
    {
        De_ = dict.lookupOrDefault<dimensionedScalar>("De", De_);
    }
    if (dict.found("gasMetalExchangeCoeff"))
    {
        if (dict.isDict("gasMetalExchangeCoeff"))
        {
            gasMetalExchangeFunction_.reset
            (
                Function1<scalar>::New
                (
                    "gasMetalExchangeCoeff",
                    dict.subDict("gasMetalExchangeCoeff")
                ).ptr()
            );
        }
        else
        {
            gasMetalExchangeCoeff_ =
                dict.lookupOrDefault<dimensionedScalar>
                (
                    "gasMetalExchangeCoeff",
                    gasMetalExchangeCoeff_
                );
        }
    }
    if (!validateParameters())
    {
        FatalErrorInFunction
            << "Invalid model parameters"
            << abort(FatalError);
    }
    // Register fields if not already present
    if (!mesh_.foundObject<volScalarField>("Te"))
    {
        mesh_.objectRegistry::store(new volScalarField(Te_));
    }
    if (!mesh_.foundObject<volScalarField>("Tl"))
    {
        mesh_.objectRegistry::store(new volScalarField(Tl_));
    }
    // Initialize energy tracking
    updateEnergyTracking();
    if (verbose)
    {
        Info<< "Two-temperature model initialized:" << nl
            << "  Ce = " << Ce_.value() << " J/m³/K" << nl
            << "  Cl = " << Cl_.value() << " J/m³/K" << nl
            << "  G = " << G_.value() << " W/m³/K" << nl
            << "  De = " << De_.value() << " m²/s" << nl
            << "  Gas-metal exchange coeff = "
            << gasMetalExchangeCoeff_.value() << " W/m³/K" << nl            
            << "  Ambient temperature = " << ambientTemperature_.value() << " K" << nl
            << "  Metal fraction field = " << metalFraction_.name() << endl;
    }
}
twoTemperatureModel::~twoTemperatureModel()
{
    if (mesh_.foundObject<volScalarField>("Te"))
    {
        mesh_.objectRegistry::checkOut("Te");
    }
    if (mesh_.foundObject<volScalarField>("Tl"))
    {
        mesh_.objectRegistry::checkOut("Tl");
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
    bool valid = true;
       // Check dimensions
    if (Te_.dimensions() != dimTemperature || 
        Tl_.dimensions() != dimTemperature)
    {
        FatalErrorInFunction
            << "Invalid temperature dimensions"
            << abort(FatalError);
        valid = false;
    }
    if (Ce_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        Cl_.dimensions() != dimEnergy/dimVolume/dimTemperature ||
        G_.dimensions() != dimEnergy/dimVolume/dimTime/dimTemperature ||
        De_.dimensions() != dimLength*dimLength/dimTime)
    {
        FatalErrorInFunction
            << "Invalid material property dimensions"
            << abort(FatalError);
        valid = false;
    }
    // Check property values
    if (Ce_.value() <= 0 || Cl_.value() <= 0 || G_.value() <= 0 || De_.value() <= 0)
    {
        FatalErrorInFunction
            << "Non-positive material properties detected"
            << abort(FatalError);
        valid = false;
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
        valid = false;
    }
    // Check temperature values
    forAll(Te_, cellI)
    {
        if (Te_[cellI] < 0 || Tl_[cellI] < 0)
        {
            FatalErrorInFunction
                << "Negative temperature detected at cell " << cellI
                << abort(FatalError);
            valid = false;
            break;
        }
    }
    return valid;
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
    const dimensionedScalar cutoffDim("metalCutoff", dimless, cutoff);
    return pos(metalFraction_ - cutoffDim);
}
tmp<volScalarField> twoTemperatureModel::clampedMetalFraction() const
{
    const scalar metalFractionFloor =
        dict_.lookupOrDefault<scalar>("metalFractionFloor", 1e-6);
    const dimensionedScalar metalFloor
    (
        "metalFractionFloor",
        dimless,
        Foam::max(metalFractionFloor, VSMALL)
    );
    const dimensionedScalar metalCeiling("metalFractionCeiling", dimless, 1.0);

    return Foam::min(Foam::max(metalFraction_, metalFloor), metalCeiling);
}

void twoTemperatureModel::applyTemperatureBounds
(
    const volScalarField& activeMask,
    const dimensionedScalar& minTemp,
    const dimensionedScalar& maxTemp,
    const dimensionedScalar& ambient
)
{
    tmp<volScalarField> tInactiveMask = scalar(1) - activeMask;
    const volScalarField& inactiveMask = tInactiveMask();

    tmp<volScalarField> tBoundTe = Foam::max(Foam::min(Te_, maxTemp), minTemp);
    tmp<volScalarField> tBoundTl = Foam::max(Foam::min(Tl_, maxTemp), minTemp);
    const volScalarField& boundTe = tBoundTe();
    const volScalarField& boundTl = tBoundTl();

    Te_ = activeMask*boundTe + inactiveMask*ambient;
    Tl_ = activeMask*boundTl + inactiveMask*ambient;

    Te_.correctBoundaryConditions();
    Tl_.correctBoundaryConditions();
}

void twoTemperatureModel::solveLatticeEquation
(
    const volScalarField& metalEff,
    const volScalarField& G,
    const volScalarField& klField,
    const volScalarField& phaseChangeRelaxCoeff,
    const volScalarField& phaseChangeSource,
    const volScalarField& TlOld,
    const volScalarField& gasMetalHeatFlux
)
{
    fvScalarMatrix TlEqn
    (
        fvm::ddt(metalEff*Cl_, Tl_)
      - fvm::laplacian(metalEff*klField, Tl_)
      + fvm::Sp(metalEff*G, Tl_)
      + fvm::Sp(metalEff*Cl_*phaseChangeRelaxCoeff, Tl_)
     ==
        metalFraction_*
        (
            G*Te_
          + Cl_*(phaseChangeSource + phaseChangeRelaxCoeff*TlOld)
        )
      + gasMetalHeatFlux
    );

    TlEqn.relax();
    TlEqn.solve(mesh_.solver("Tl"));
    Tl_.correctBoundaryConditions();
}

void twoTemperatureModel::solveElectronEquation
(
    const volScalarField& metalEff,
    const volScalarField& Ce,
    const volScalarField& ke,
    const volScalarField& G,
    const volScalarField& laserSource
)
{
    fvScalarMatrix TeEqn
    (
        fvm::ddt(metalEff*Ce, Te_)
      - fvm::laplacian(metalEff*ke, Te_)
      + fvm::Sp(metalEff*G, Te_)
     ==
        metalFraction_*(laserSource + G*Tl_)
    );

    TeEqn.relax();
    TeEqn.solve(mesh_.solver("Te"));
    Te_.correctBoundaryConditions();
}

dimensionedScalar twoTemperatureModel::couplingEnergy
(
    const volScalarField& gasMetalHeatFluxMasked,
    const dimensionedScalar& dtDim
) const
{
    return fvc::domainIntegrate(gasMetalHeatFluxMasked)*dtDim;
}

dimensionedScalar twoTemperatureModel::currentTotalEnergy() const
{
    tmp<volScalarField> tCe = electronHeatCapacity();
    const volScalarField& CeField = tCe();
   tmp<volScalarField> tMetal = clampedMetalFraction();
    const volScalarField& metalEff = tMetal();
    return fvc::domainIntegrate(metalEff*(CeField*Te_ + Cl_*Tl_));
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
        fvc::domainIntegrate(metalEff*(CeAfter*Te_));
    const dimensionedScalar latticeEnergyAfter =
        fvc::domainIntegrate(metalEff*(Cl_*Tl_));

    static dimensionedScalar cumulativeLaserEnergy
    (
        "cumulativeLaserEnergy",        
        dimEnergy,
        0.0
    );
    cumulativeLaserEnergy += laserEnergy;

    Info<< "Energy diagnostics:" << nl
        << "  Electron energy before: "
        << electronEnergyBefore.value() << " J" << nl
        << "  Lattice energy before: "
        << latticeEnergyBefore.value() << " J" << nl
        << "  Laser energy input: " << laserEnergy.value() << " J" << nl
        << "  Phase-change energy input: "
        << phaseChangeEnergy.value() << " J" << nl
        << "  Gas-metal coupling energy input: "
        << couplingEnergy.value() << " J" << nl
        << "  Electron energy after: "
        << electronEnergyAfter.value() << " J" << nl
        << "  Lattice energy after: "
        << latticeEnergyAfter.value() << " J" << nl
        << "  Total metal energy: "
        << (electronEnergyAfter + latticeEnergyAfter).value() << " J" << nl
        << "  Cumulative laser energy: "
        << cumulativeLaserEnergy.value() << " J" << endl;
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
        << "Source energy: " << totalEnergyInput.value() << " J" << nl
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
    lastTotalEnergy_ = currentEnergy;
    energyInitialized_ = true;
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

    const dimensionedScalar minTemp
    (
        "minTemp",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("minTe", 300.0)
    );
    const dimensionedScalar maxTemp
    (
        "maxTemp",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("maxTe", 3500.0)
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

    tmp<volScalarField> tActiveMask = metalActiveMask(metalCutoff);
    const volScalarField& activeMask = tActiveMask();
    applyTemperatureBounds(activeMask, minTemp, maxTemp, ambientDim);

    dimensionedScalar previousEnergy("previousEnergy", dimEnergy, 0.0);
    if (energyInitialized_)
    {
        previousEnergy = lastTotalEnergy_;
    }
    else
    {
        previousEnergy = currentTotalEnergy();
    }

    tmp<volScalarField> tCouplingMask =
        metalActiveMask(Foam::max(metalFractionFloor, VSMALL));
    const volScalarField& couplingMask = tCouplingMask();
    tmp<volScalarField> tGasMetalHeatFluxMasked = couplingMask*gasMetalHeatFlux;
    const volScalarField& gasMetalHeatFluxMasked = tGasMetalHeatFluxMasked();

    const dimensionedScalar dtDim = mesh_.time().deltaT();
    const dimensionedScalar laserEnergy =
        fvc::domainIntegrate(metalEff*laserSource)*dtDim;
    const dimensionedScalar phaseChangeEnergy =
        fvc::domainIntegrate(metalEff*(Cl_*phaseChangeSource))*dtDim;
    const dimensionedScalar couplingEnergyValue =
        couplingEnergy(gasMetalHeatFluxMasked, dtDim);
    const dimensionedScalar totalEnergyInput =
        laserEnergy + phaseChangeEnergy + couplingEnergyValue;

    dimensionedScalar electronEnergyBefore("electronEnergyBefore", dimEnergy, 0.0);
    dimensionedScalar latticeEnergyBefore("latticeEnergyBefore", dimEnergy, 0.0);

    if (energyDiagnostics)
    {
        tmp<volScalarField> tCeInitial = electronHeatCapacity();
        const volScalarField& CeInitial = tCeInitial();
        electronEnergyBefore =
            fvc::domainIntegrate(metalEff*(CeInitial*Te_));
        latticeEnergyBefore =
            fvc::domainIntegrate(metalEff*(Cl_*Tl_));
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
        tmp<volScalarField> tG = electronPhononCoupling();
        const volScalarField& G = tG();

        tmp<volScalarField> tkl = kl();
        const volScalarField& klField = tkl();

        solveLatticeEquation
        (
            metalEff,
            G,
            klField,
            phaseChangeRelaxCoeff,
            phaseChangeSource,
            TlOld,
            gasMetalHeatFluxMasked
        );

        tmp<volScalarField> tke = electronThermalConductivity();
        const volScalarField& keField = tke();

        tmp<volScalarField> tCe = electronHeatCapacity();
        const volScalarField& CeField = tCe();

        solveElectronEquation
        (
            metalEff,
            CeField,
            keField,
            G,
            laserSource
        );

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

    applyTemperatureBounds(activeMask, minTemp, maxTemp, ambientDim);
    residual = gMax(mag(Te_ - Tl_)().internalField());

    tmp<volScalarField> tCeFinal = electronHeatCapacity();
    const volScalarField& CeFinal = tCeFinal();
    const dimensionedScalar currentEnergy =
        fvc::domainIntegrate(metalEff*(CeFinal*Te_ + Cl_*Tl_));

    scalar energyError = 0.0;
    if
    (
        !checkEnergyConservation
        (
            previousEnergy,
            totalEnergyInput,
            currentEnergy,
            energyError
        )
    )
    {
        reportEnergyViolation
        (
            previousEnergy,
            totalEnergyInput,
            currentEnergy,
            energyError,
            dtDim.value(),
            energyAudit
        );
    }

    updateEnergyTracking(currentEnergy);

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
    // Create temporary field for electron thermal conductivity
    tmp<volScalarField> tke
    (
        new volScalarField
        (
            IOobject
            (
                "ke",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("ke", dimPower/dimLength/dimTemperature, 1.0)
        )
    );
       // Get reference to field for modification
    volScalarField& ke = tke.ref(); // Calculate conductivity using cell-wise electron heat capacity
    // ke = CeField * De_, allowing spatially varying Ce
    tmp<volScalarField> tCe = electronHeatCapacity();
    const volScalarField& CeField = tCe();
    forAll(ke, cellI)
    {
        ke[cellI] = (CeField[cellI] * De_).value();
    }
    return tke;
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
                IOobject::NO_WRITE
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
            CeField[cellI] = CeFunction_->value(Te_[cellI]);
        }
    }
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
                IOobject::NO_WRITE
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
            GField[cellI] = GFunction_->value(Te_[cellI]);
        }
    }
    return tG;
}
tmp<volScalarField> twoTemperatureModel::gasMetalExchangeCoeffField() const
{
    tmp<volScalarField> tCoeff
    (
        new volScalarField
        (
            IOobject("gasMetalExchangeCoeffField", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_,
            gasMetalExchangeCoeff_
        )
    );
    volScalarField& coeff = tCoeff.ref();
    if (gasMetalExchangeFunction_.valid())
    {
        forAll(coeff, cellI)
        {
            coeff[cellI] = gasMetalExchangeFunction_->value(Tl_[cellI]);
        }
    }
    coeff *= metalFraction_;
    coeff *= (scalar(1) - metalFraction_);
    const scalar maxExchange =
        dict_.lookupOrDefault<scalar>("maxGasMetalExchangeCoeff", 1e18);
    const scalar minExchange =
        dict_.lookupOrDefault<scalar>("minGasMetalExchangeCoeff", 0.0);
    const dimensionedScalar minExchangeDim
    (
        "minGasMetalExchangeCoeff",
        coeff.dimensions(),
        minExchange
    );
    const dimensionedScalar maxExchangeDim
    (
        "maxGasMetalExchangeCoeff",
        coeff.dimensions(),
        maxExchange
    );
    coeff = Foam::max(Foam::min(coeff, maxExchangeDim), minExchangeDim);
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
    if (verbose)
    {
        Info<< "Two-temperature model:" << nl
            << "Parameters:" << nl
            << "  Ce = " << Ce_.value() << " J/m³/K" << nl
            << "  Cl = " << Cl_.value() << " J/m³/K" << nl
            << "  G = " << G_.value() << " W/m³/K" << nl
            << "Field statistics:" << nl
            << "  Te range: " << min(Te_).value() << " - " << max(Te_).value() << " K" << nl
            << "  Tl range: " << min(Tl_).value() << " - " << max(Tl_).value() << " K" << nl
            << "  Mean Te: " << average(Te_).value() << " K" << nl
            << "  Mean Tl: " << average(Tl_).value() << " K" << nl;
    }
    if (energyInitialized_)
    {
        tmp<volScalarField> tCeW = electronHeatCapacity();
        const volScalarField& CeW = tCeW();
        tmp<volScalarField> tMetal = clampedMetalFraction();
        const volScalarField& metalEff = tMetal();        
        dimensionedScalar currentEnergy =
            fvc::domainIntegrate(metalEff*(CeW*Te_ + Cl_*Tl_));
        scalar energyError = mag
        (
            (currentEnergy.value() - lastTotalEnergy_.value())/
            (mag(lastTotalEnergy_.value()) + SMALL)
        );
        if (verbose)
        {
            Info<< "Energy conservation:" << nl
                << "  Current total energy: " << currentEnergy.value() << " J" << nl
                << "  Energy error: " << energyError * 100 << " %" << endl;
        }
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
            IOobject("kl", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
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
