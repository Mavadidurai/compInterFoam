/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "plasmaIonizationModel.H"
#include "Pstream.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include <cmath>

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plasmaIonizationModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plasmaIonizationModel::plasmaIonizationModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    ionizationEnergy_
    (
        "ionizationEnergy",
        dimensionSet(1, 2, -2, 0, 0),  // [J]
        dict.lookupOrDefault<scalar>("ionizationEnergy", 6.82 * 1.602e-19)  // Ti: 6.82 eV
    ),
    T_ionization_
    (
        "T_ionization",
        dimTemperature,
        dict.lookupOrDefault<scalar>("T_ionization", 30000.0)  // ~30,000 K threshold
    ),
    k_B_
    (
        "k_B",
        dimensionSet(1, 2, -2, -1, 0),
        physicoChemical::k.value()  // 1.38064852e-23 J/K
    ),
    m_e_
    (
        "m_e",
        dimensionSet(1, 0, 0, 0, 0),
        electromagnetic::me.value()  // 9.10938356e-31 kg
    ),
    m_atom_
    (
        "m_atom",
        dimensionSet(1, 0, 0, 0, 0),
        dict.lookupOrDefault<scalar>("atomicMass", 7.95e-26)  // Ti: 47.867 amu
    ),
    n_atom_
    (
        "n_atom",
        dimensionSet(0, -3, 0, 0, 0),
        dict.lookupOrDefault<scalar>("atomicNumberDensity", 5.68e28)  // Ti: ~5.68e28 atoms/m^3
    ),
    enabled_
    (
        dict.lookupOrDefault<Switch>("enablePlasmaModel", true)
    ),
    ionizationDegree_
    (
        IOobject
        (
            "ionizationDegree",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    n_electron_
    (
        IOobject
        (
            "n_electron",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -3, 0, 0, 0), 0.0)
    ),
    plasmaPressure_
    (
        IOobject
        (
            "plasmaPressure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
    plasmaShielding_
    (
        IOobject
        (
            "plasmaShielding",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    )
{
    if (enabled_ && Pstream::master())
    {
        Info<< "Plasma ionization model enabled:" << nl
            << "    Ionization energy = " << ionizationEnergy_.value()/1.602e-19 << " eV" << nl
            << "    T_ionization = " << T_ionization_.value() << " K" << nl
            << "    Atomic mass = " << m_atom_.value() << " kg" << nl
            << "    Atomic density = " << n_atom_.value() << " atoms/m^3" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::plasmaIonizationModel::~plasmaIonizationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::plasmaIonizationModel::sahaIonization
(
    const scalar T,
    const scalar n_atoms
) const
{
    if (T < T_ionization_.value())
    {
        return 0.0;
    }

    const scalar k_B = k_B_.value();
    const scalar E_ion = ionizationEnergy_.value();
    const scalar m_e = m_e_.value();

    // Saha equation (simplified, single ionization):
    // α^2/(1-α) = (2πm_e k_B T/h^2)^(3/2) * 2/n * exp(-E_ion/k_B T)
    //
    // For T >> T_ionization, α → 1
    // For T ~ T_ionization, α ~ 0.01 - 0.1

    const scalar h = 6.62607015e-34;  // Planck constant [J·s]
    const scalar pi = mathematical::pi;

    // Thermal de Broglie wavelength squared
    const scalar lambda_dB_sq = sqr(h) / (2.0 * pi * m_e * k_B * T);

    // Partition function ratio (assume g_ion/g_atom ~ 2)
    const scalar g_ratio = 2.0;

    // Saha coefficient
    const scalar K_saha = g_ratio * pow(2.0 / (lambda_dB_sq * n_atoms), 1.5)
                        * exp(-E_ion / (k_B * T));

    // Solve quadratic: α^2 + K_saha * α - K_saha = 0
    // α = (-K_saha + sqrt(K_saha^2 + 4*K_saha)) / 2

    const scalar discriminant = sqr(K_saha) + 4.0 * K_saha;
    const scalar alpha = max
    (
        0.0,
        min
        (
            1.0,
            (-K_saha + sqrt(discriminant)) / 2.0
        )
    );

    return alpha;
}


void Foam::plasmaIonizationModel::update
(
    const volScalarField& T,
    const volScalarField& rho,
    const volScalarField& alpha1
)
{
    if (!enabled_)
    {
        return;
    }

    // Reset fields
    ionizationDegree_.primitiveFieldRef() = 0.0;
    n_electron_.primitiveFieldRef() = 0.0;
    plasmaPressure_.primitiveFieldRef() = 0.0;
    plasmaShielding_.primitiveFieldRef() = 0.0;

    const scalar T_ion = T_ionization_.value();
    const scalar k_B = k_B_.value();
    const scalar n_atom = n_atom_.value();
    const scalar m_a = m_atom_.value();

    const scalarField& TField = T.primitiveField();
    const scalarField& rhoField = rho.primitiveField();
    const scalarField& alpha1Field = alpha1.primitiveField();

    scalarField& alphaField = ionizationDegree_.primitiveFieldRef();
    scalarField& n_e_Field = n_electron_.primitiveFieldRef();
    scalarField& p_plasma_Field = plasmaPressure_.primitiveFieldRef();
    scalarField& shieldField = plasmaShielding_.primitiveFieldRef();

    label nPlasmaCells = 0;

    forAll(TField, cellI)
    {
        const scalar T_local = TField[cellI];
        const scalar rho_local = rhoField[cellI];
        const scalar alpha = alpha1Field[cellI];

        // Only apply to metal vapor/plasma regions (alpha > 0.1, T > T_ion)
        if (T_local > T_ion && alpha > 0.1)
        {
            // Estimate atom density from material density
            const scalar n_local = max(rho_local / m_a, SMALL);

            // Calculate ionization degree using Saha equation
            const scalar alpha_ion = sahaIonization(T_local, n_local);

            alphaField[cellI] = alpha_ion;

            // Electron density
            const scalar n_e = alpha_ion * n_local;
            n_e_Field[cellI] = n_e;

            // Plasma pressure: p = (n_e + n_ion) * k_B * T
            // For single ionization: n_ion = n_e
            // p = 2 * n_e * k_B * T
            const scalar p_plasma = 2.0 * n_e * k_B * T_local;
            p_plasma_Field[cellI] = p_plasma;

            // Plasma shielding: reduces laser absorption
            // Shield coefficient ~ 1 - exp(-n_e / n_critical)
            // n_critical ~ ε_0 * m_e * ω_laser^2 / e^2
            // For 343 nm (Ti:Sapphire 3rd harmonic): ω = 5.5e15 rad/s
            // n_crit ~ 1.1e27 m^-3
            const scalar n_crit = 1.1e27;
            const scalar shield = 1.0 - exp(-n_e / n_crit);
            shieldField[cellI] = shield;

            nPlasmaCells++;
        }
    }

    // Correct boundary conditions
    ionizationDegree_.correctBoundaryConditions();
    n_electron_.correctBoundaryConditions();
    plasmaPressure_.correctBoundaryConditions();
    plasmaShielding_.correctBoundaryConditions();

    // Report statistics
    if (nPlasmaCells > 0 && Pstream::master())
    {
        const scalar maxIonization = gMax(ionizationDegree_);
        const scalar maxPlasmaP = gMax(plasmaPressure_);
        const scalar maxShielding = gMax(plasmaShielding_);
        const scalar maxNe = gMax(n_electron_);

        if (mesh_.time().timeIndex() % 100 == 0)
        {
            Info<< "Plasma formation active:" << nl
                << "    Plasma cells = " << returnReduce(nPlasmaCells, sumOp<label>()) << nl
                << "    Max ionization degree = " << maxIonization << nl
                << "    Max electron density = " << maxNe << " m^-3" << nl
                << "    Max plasma pressure = " << maxPlasmaP/1e6 << " MPa" << nl
                << "    Max laser shielding = " << maxShielding*100.0 << " %" << endl;
        }
    }
}


// ************************************************************************* //
