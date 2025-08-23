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
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
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
        dict.lookupOrDefault<dimensionedScalar>
        (
            "Ce",
            dimensionedScalar
            (
                "Ce_default",
                dimEnergy/dimVolume/dimTemperature,
                210.0
            )
        ).value()
    ),
    Cl_
    (
        "Cl",
        dimEnergy/dimVolume/dimTemperature,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "Cl",
            dimensionedScalar
            (
                "Cl_default",
                dimEnergy/dimVolume/dimTemperature,
                2.3e6
            )
        ).value()
    ),
    G_
    (
        "G",
        dimEnergy/dimVolume/dimTime/dimTemperature,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "G",
            dimensionedScalar
            (
                "G_default",
                dimEnergy/dimVolume/dimTime/dimTemperature,
                5e17
            )
        ).value()
    ),
        De_
    (
        "De",
        dimLength*dimLength/dimTime,
        dict.lookupOrDefault<dimensionedScalar>
        (
            "De",
            dimensionedScalar
            (
                "De_default",
                dimLength*dimLength/dimTime,
                1e-4
            )
        ).value()
    ),
    lastTotalEnergy_
    (
        "lastTotalEnergy",
        dimEnergy,
        0.0
    ),
    energyInitialized_(false)
{
        const bool verbose =
        mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);
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
            << "  De = " << De_.value() << " m²/s" << endl;
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

bool twoTemperatureModel::validateFields() const
{
    bool valid = true;

    forAll(Te_, cellI)
    {
        if (Te_[cellI] < 0 || !std::isfinite(Te_[cellI]) ||
            Tl_[cellI] < 0 || !std::isfinite(Tl_[cellI]))
        {
            valid = false;
            break;
        }
    }

    return valid;
}

bool twoTemperatureModel::checkEnergyConservation() const
{
    if (!energyInitialized_)
    {
        return true;
    }

    dimensionedScalar currentEnergy = 
        fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_);

    scalar energyError = mag
    (
        (currentEnergy.value() - lastTotalEnergy_.value())/
        (mag(lastTotalEnergy_.value()) + SMALL)
    );

    return energyError < dict_.get<scalar>("energyTolerance");
}

void twoTemperatureModel::updateEnergyTracking() const
{
    lastTotalEnergy_ = fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_);
    energyInitialized_ = true;
}

void twoTemperatureModel::solve
(
    const volScalarField& laserSource,
    const volScalarField& phaseChangeSource
)
{
    if (!validateFields())
    {
        FatalErrorInFunction
            << "Invalid field values before solve"
            << abort(FatalError);
    }

    // Store initial energy
    updateEnergyTracking();
    if (verbose)
    {
        Info<< "max(laserSource) = " << max(laserSource).value()
            << ", max(Tl_) = " << max(Tl_).value() << endl;
    }


    // Calculate temperature-dependent properties
    volScalarField ke = electronThermalConductivity();
    volScalarField Ce = electronHeatCapacity();
    volScalarField G = electronPhononCoupling();

    // Apply strong under-relaxation for stability
    scalar relaxFactor = 0.3;

    // First solve lattice temperature - more stable to do this first
    // Create a more stable matrix system for the lattice temperature
    fvScalarMatrix TlEqn
    (
        fvm::ddt(Cl_, Tl_)
      - fvm::laplacian(kl(), Tl_)
      + fvm::Sp(G, Tl_)
     ==
        G*Te_
      + Cl_*phaseChangeSource

    );

    // Under-relax the lattice equation for stability
    TlEqn.relax(relaxFactor);
    
    // Use more robust solver settings for lattice temperature
    dictionary latticeDict;
    latticeDict.add("solver", "PBiCGStab");
    latticeDict.add("preconditioner", "DIC");
    latticeDict.add("tolerance", 1e-8);
    latticeDict.add("relTol", 0.01);
    latticeDict.add("maxIter", 500);
    
    // Perform lattice temperature solve; any additional iterations are
    // handled by the solver controls (maxIter, tolerance, etc.)
    TlEqn.solve(latticeDict);

    // Create a more stable matrix system for the electron temperature
    fvScalarMatrix TeEqn
    (
        fvm::ddt(Ce, Te_)
      - fvm::laplacian(ke, Te_)
      + fvm::Sp(G, Te_)
     ==
        laserSource
      + G*Tl_
    );


    // Apply stronger under-relaxation to electron equation
    TeEqn.relax(relaxFactor);
    
    // Use customized robust solver settings
    dictionary electronDict;
    electronDict.add("solver", "PBiCGStab");
    electronDict.add("preconditioner", "DIC");
    electronDict.add("tolerance", 1e-6);
    electronDict.add("relTol", 0.01);  // Looser relative tolerance
    electronDict.add("maxIter", 1000);
    
    TeEqn.solve(electronDict);

    // Apply temperature bounds from dictionary
    dimensionedScalar minTemp
    (
        "minTemp",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("minTe", 300.0)
    );
    dimensionedScalar maxTemp
    (
        "maxTemp",
        dimTemperature,
        dict_.lookupOrDefault<scalar>("maxTe", 3500.0)
    );

    // Bound fields safely
    forAll(Te_, cellI)
    {
        Te_[cellI] = max(min(Te_[cellI], maxTemp.value()), minTemp.value());
        Tl_[cellI] = max(min(Tl_[cellI], maxTemp.value()), minTemp.value());
    }

    Te_.correctBoundaryConditions();
    Tl_.correctBoundaryConditions();
    
    // Check energy conservation
    if (!checkEnergyConservation())
    {
        dimensionedScalar currentEnergy = fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_);
        scalar energyError = mag
        (
            (currentEnergy.value() - lastTotalEnergy_.value())/
            (mag(lastTotalEnergy_.value()) + SMALL)
        );
        
        // Characteristic time-scales based on Ce/G and Cl/G
        scalar dt = mesh_.time().deltaTValue();
        scalar tauE = Ce_.value()/G_.value();
        scalar tauL = Cl_.value()/G_.value();

        WarningInFunction
            << "Energy conservation violation detected" << nl
            << "Error = " << energyError * 100 << " %" << nl
            << "Previous energy: " << lastTotalEnergy_.value() << " J" << nl
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

    // Update energy tracking
    updateEnergyTracking();

    // Diagnostics: track cumulative laser energy versus lattice/electron energy
    {
        static dimensionedScalar cumulativeLaserEnergy
        (
            "cumulativeLaserEnergy",
            dimEnergy,
            0.0
        );

        dimensionedScalar laserEnergyThisStep =
            fvc::domainIntegrate(laserSource) * mesh_.time().deltaT();
        cumulativeLaserEnergy += laserEnergyThisStep;

        dimensionedScalar electronEnergy = fvc::domainIntegrate(Ce_*Te_);
        dimensionedScalar latticeEnergy  = fvc::domainIntegrate(Cl_*Tl_);

        Info<< "Energy diagnostics:" << nl
            << "  Cumulative laser energy: "
            << cumulativeLaserEnergy.value() << " J" << nl
            << "  Electron energy: " << electronEnergy.value() << " J" << nl
            << "  Lattice energy: " << latticeEnergy.value() << " J" << nl
            << "  Total energy: "
            << (electronEnergy + latticeEnergy).value() << " J" << endl;
    }

    // Report solution statistics with more detail
    volScalarField tempDiff = mag(Te_ - Tl_);

    if (verbose)
    {
        Info<< "Two-temperature solve:" << nl
            << "  Te range: " << min(Te_).value() << " - " << max(Te_).value() << " K" << nl
            << "  Tl range: " << min(Tl_).value() << " - " << max(Tl_).value() << " K" << nl
            << "  Max temperature difference: " << max(tempDiff).value() << " K" << nl
            << "  Mean Te: " << gAverage(Te_) << " K" << nl
            << "  Mean Tl: " << gAverage(Tl_) << " K" << endl;
    }
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
    volScalarField& ke = tke.ref();
    // Calculate conductivity using cell-wise electron heat capacity
    // ke = CeField * De_, allowing spatially varying Ce
    const volScalarField CeField = electronHeatCapacity();
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
    );

    return tG;
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
        dimensionedScalar currentEnergy = fvc::domainIntegrate(Ce_*Te_ + Cl_*Tl_);
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
    // Electronic thermal conductivity derived from constant diffusivity
    return electronThermalConductivity();
}
tmp<volScalarField> Foam::twoTemperatureModel::kl() const
{
    // Return lattice thermal conductivity
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
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("kl", dimPower/dimLength/dimTemperature, 0.0)
        )
    );

    volScalarField& kl = tkl.ref();

    // Calculate temperature-dependent lattice thermal conductivity
    // Using Drude model with phonon contribution
    forAll(mesh_.C(), cellI)
    {
        scalar Tl = Tl_[cellI];
        kl[cellI] = Cl_.value() * sqrt(Tl) / (3.0 * G_.value());
        
        // Apply temperature dependent correction
        if (Tl > 1000.0)
        {
            kl[cellI] *= pow(1000.0/Tl, 0.5);  // High temperature correction
        }
    }

    return tkl;
}
} // End namespace Foam
