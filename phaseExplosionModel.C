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

#include "phaseExplosionModel.H"
#include "Pstream.H"
#include "mathematicalConstants.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseExplosionModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseExplosionModel::phaseExplosionModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    T_critical_
    (
        "T_critical",
        dimTemperature,
        dict.lookupOrDefault<scalar>("T_critical", 9500.0)  // Ti critical temp
    ),
    T_spinodal_
    (
        "T_spinodal",
        dimTemperature,
        dict.lookupOrDefault<scalar>("T_spinodal", 0.9*T_critical_.value())
    ),
    p_critical_
    (
        "p_critical",
        dimPressure,
        dict.lookupOrDefault<scalar>("p_critical", 1.2e8)  // Ti: ~120 MPa
    ),
    tau_explosion_
    (
        "tau_explosion",
        dimTime,
        dict.lookupOrDefault<scalar>("tau_explosion", 1e-12)  // ps timescale
    ),
    explosionMultiplier_
    (
        dict.lookupOrDefault<scalar>("explosionMultiplier", 100.0)
    ),
    enabled_
    (
        dict.lookupOrDefault<Switch>("enablePhaseExplosion", true)
    ),
    explosiveMassSource_
    (
        IOobject
        (
            "explosiveMassSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)
    ),
    explosivePressure_
    (
        IOobject
        (
            "explosivePressure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
    explosionIndicator_
    (
        IOobject
        (
            "explosionIndicator",
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
        Info<< "Phase explosion model enabled:" << nl
            << "    T_critical   = " << T_critical_.value() << " K" << nl
            << "    T_spinodal   = " << T_spinodal_.value() << " K" << nl
            << "    p_critical   = " << p_critical_.value()/1e6 << " MPa" << nl
            << "    tau_explosion = " << tau_explosion_.value() << " s" << nl
            << "    explosionMultiplier = " << explosionMultiplier_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseExplosionModel::~phaseExplosionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseExplosionModel::update
(
    const volScalarField& T,
    const volScalarField& alpha1,
    const volScalarField& rho
)
{
    if (!enabled_)
    {
        return;
    }

    // Reset fields
    explosiveMassSource_.primitiveFieldRef() = 0.0;
    explosivePressure_.primitiveFieldRef() = 0.0;
    explosionIndicator_.primitiveFieldRef() = 0.0;

    const scalar T_spin = T_spinodal_.value();
    const scalar T_crit = T_critical_.value();
    const scalar tau = tau_explosion_.value();

    const scalarField& TField = T.primitiveField();
    const scalarField& alpha1Field = alpha1.primitiveField();
    const scalarField& rhoField = rho.primitiveField();

    scalarField& massSourceField = explosiveMassSource_.primitiveFieldRef();
    scalarField& pressureField = explosivePressure_.primitiveFieldRef();
    scalarField& indicatorField = explosionIndicator_.primitiveFieldRef();

    label nExplosiveCells = 0;

    forAll(TField, cellI)
    {
        const scalar T_local = TField[cellI];
        const scalar alpha = alpha1Field[cellI];
        const scalar rho_local = rhoField[cellI];

        // Only apply to liquid metal (alpha > 0.5) above spinodal
        if (T_local > T_spin && alpha > 0.5)
        {
            // Normalized superheat: xi = (T - T_spin) / (T_crit - T_spin)
            const scalar xi = min
            (
                (T_local - T_spin) / max(T_crit - T_spin, SMALL),
                1.0
            );

            // Explosion indicator (0 to 1)
            indicatorField[cellI] = xi;

            // Enhanced mass source: exponentially increasing with superheat
            // J_explosion ~ rho / tau * exp(10*xi)
            const scalar explosionRate = (rho_local / tau) * exp(10.0 * xi);
            massSourceField[cellI] = explosionRate;

            // Explosive pressure: scales with superheat
            // p_explosion ~ p_critical * xi^2
            pressureField[cellI] = p_critical_.value() * sqr(xi);

            nExplosiveCells++;
        }
    }

    // Correct boundary conditions
    explosiveMassSource_.correctBoundaryConditions();
    explosivePressure_.correctBoundaryConditions();
    explosionIndicator_.correctBoundaryConditions();

    // Report statistics
    if (nExplosiveCells > 0 && Pstream::master())
    {
        const scalar maxT = gMax(T);
        const scalar maxExplosionP = gMax(explosivePressure_);
        const scalar maxMassSource = gMax(explosiveMassSource_);

        if (mesh_.time().timeIndex() % 100 == 0)
        {
            Info<< "Phase explosion active:" << nl
                << "    Explosive cells = " << returnReduce(nExplosiveCells, sumOp<label>()) << nl
                << "    Max T = " << maxT << " K" << nl
                << "    Max explosion pressure = " << maxExplosionP/1e6 << " MPa" << nl
                << "    Max mass source = " << maxMassSource << " kg/m^3/s" << endl;
        }
    }
}


Foam::scalar Foam::phaseExplosionModel::explosionFluxMultiplier
(
    const scalar T
) const
{
    if (!enabled_ || T < T_spinodal_.value())
    {
        return 1.0;
    }

    const scalar T_spin = T_spinodal_.value();
    const scalar T_crit = T_critical_.value();

    // Normalized superheat
    const scalar xi = min
    (
        (T - T_spin) / max(T_crit - T_spin, SMALL),
        1.0
    );

    // Exponential enhancement above spinodal
    // multiplier = 1 + (explosionMultiplier - 1) * xi^2
    return 1.0 + (explosionMultiplier_ - 1.0) * sqr(xi);
}


Foam::scalar Foam::phaseExplosionModel::explosionPressure
(
    const scalar T,
    const scalar rho
) const
{
    if (!enabled_ || T < T_spinodal_.value())
    {
        return 0.0;
    }

    const scalar T_spin = T_spinodal_.value();
    const scalar T_crit = T_critical_.value();

    const scalar xi = min
    (
        (T - T_spin) / max(T_crit - T_spin, SMALL),
        1.0
    );

    return p_critical_.value() * sqr(xi);
}


// ************************************************************************* //
