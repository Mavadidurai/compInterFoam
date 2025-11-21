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

#include "enhancedLIFTPhysics.H"
#include "twoPhaseMixtureThermo.H"
#include "fvc.H"
#include "Pstream.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"
#include "electromagneticConstants.H"
#include "physicoChemicalConstants.H"
#include "atomicConstants.H"
#include <cmath>

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(enhancedLIFTPhysics, 0);
}


// * * * * * * * * * * * PhaseExplosionData Constructors * * * * * * * * * * //

Foam::enhancedLIFTPhysics::PhaseExplosionData::PhaseExplosionData
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    T_critical_
    (
        "T_critical",
        dimTemperature,
        dict.lookupOrDefault<scalar>("T_critical", 9500.0)
    ),
    T_spinodal_
    (
        "T_spinodal",
        dimTemperature,
        dict.lookupOrDefault<scalar>("T_spinodal", 0.9*dict.lookupOrDefault<scalar>("T_critical", 9500.0))
    ),
    p_critical_
    (
        "p_critical",
        dimPressure,
        dict.lookupOrDefault<scalar>("p_critical", 1.2e8)
    ),
    tau_explosion_
    (
        "tau_explosion",
        dimTime,
        dict.lookupOrDefault<scalar>("tau_explosion", 1e-12)
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
        dimensionedScalar
        (
            "zero",
            dimensionSet(1, -3, -1, 0, 0, 0, 0),
            0.0
        )
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


// * * * * * * * * * * * PlasmaData Constructors  * * * * * * * * * * * * * //

Foam::enhancedLIFTPhysics::PlasmaData::PlasmaData
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ionizationEnergy_
    (
        "ionizationEnergy",
        dimensionSet(1, 2, -2, 0, 0, 0, 0),
        dict.lookupOrDefault<scalar>("ionizationEnergy", 6.82 * 1.602e-19)
    ),
    T_ionization_
    (
        "T_ionization",
        dimTemperature,
        dict.lookupOrDefault<scalar>("T_ionization", 30000.0)
    ),
    k_B_
    (
        "k_B",
        dimensionSet(1, 2, -2, -1, 0, 0, 0),
        1.38064852e-23  // Boltzmann constant [J/K]
    ),
    m_e_
    (
        "m_e",
        dimensionSet(1, 0, 0, 0, 0, 0, 0),
        9.10938356e-31  // Electron mass [kg]
    ),
    m_atom_
    (
        "m_atom",
        dimensionSet(1, 0, 0, 0, 0, 0, 0),
        dict.lookupOrDefault<scalar>("atomicMass", 7.95e-26)
    ),
    n_atom_
    (
        "n_atom",
        dimensionSet(0, -3, 0, 0, 0, 0, 0),
        dict.lookupOrDefault<scalar>("atomicNumberDensity", 5.68e28)
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
        dimensionedScalar
        (
            "zero",
            dimensionSet(0, -3, 0, 0, 0, 0, 0),
            0.0
        )
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


// * * * * * * * * * * * BreakupData Constructors * * * * * * * * * * * * * //

Foam::enhancedLIFTPhysics::BreakupData::BreakupData
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    We_critical_
    (
        dict.lookupOrDefault<scalar>("We_critical", 10.0)
    ),
    minJetDiameter_
    (
        "minJetDiameter",
        dimLength,
        dict.lookupOrDefault<scalar>("minJetDiameter", 1e-7)
    ),
    maxJetDiameter_
    (
        "maxJetDiameter",
        dimLength,
        dict.lookupOrDefault<scalar>("maxJetDiameter", 10e-6)
    ),
    breakupTimeCoeff_
    (
        dict.lookupOrDefault<scalar>("breakupTimeCoeff", 1.0)
    ),
    enabled_
    (
        dict.lookupOrDefault<Switch>("enableDropletBreakup", true)
    ),
    WeberNumber_
    (
        IOobject
        (
            "WeberNumber",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    breakupIndicator_
    (
        IOobject
        (
            "breakupIndicator",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    dropletDiameter_
    (
        IOobject
        (
            "dropletDiameter",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    breakupRate_
    (
        IOobject
        (
            "breakupRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    )
{
    if (enabled_ && Pstream::master())
    {
        Info<< "Droplet breakup model enabled:" << nl
            << "    We_critical = " << We_critical_ << nl
            << "    Min jet diameter = " << minJetDiameter_.value()*1e6 << " μm" << nl
            << "    Max jet diameter = " << maxJetDiameter_.value()*1e6 << " μm" << nl
            << "    Breakup time coefficient = " << breakupTimeCoeff_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enhancedLIFTPhysics::enhancedLIFTPhysics
(
    const fvMesh& mesh,
    const twoPhaseMixtureThermo& mixture,
    const dictionary& dict
)
:
    mesh_(mesh),
    mixture_(mixture),
    phaseExplosion_(nullptr),
    plasma_(nullptr),
    breakup_(nullptr)
{
    // Initialize phase explosion model
    if (dict.found("phaseExplosionCoeffs"))
    {
        phaseExplosion_.reset
        (
            new PhaseExplosionData
            (
                mesh,
                dict.subDict("phaseExplosionCoeffs")
            )
        );
    }

    // Initialize plasma ionization model
    if (dict.found("plasmaIonizationCoeffs"))
    {
        plasma_.reset
        (
            new PlasmaData
            (
                mesh,
                dict.subDict("plasmaIonizationCoeffs")
            )
        );
    }

    // Initialize droplet breakup model
    if (dict.found("dropletBreakupCoeffs"))
    {
        breakup_.reset
        (
            new BreakupData
            (
                mesh,
                dict.subDict("dropletBreakupCoeffs")
            )
        );
    }

    if (Pstream::master())
    {
        Info<< "Enhanced LIFT physics initialized:" << nl
            << "    Phase explosion: " << (phaseExplosionEnabled() ? "ON" : "OFF") << nl
            << "    Plasma ionization: " << (plasmaEnabled() ? "ON" : "OFF") << nl
            << "    Droplet breakup: " << (breakupEnabled() ? "ON" : "OFF") << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::enhancedLIFTPhysics::~enhancedLIFTPhysics()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::enhancedLIFTPhysics::sahaIonization
(
    const scalar T,
    const scalar n_atoms
) const
{
    if (!plasma_.valid() || T < plasma_->T_ionization_.value())
    {
        return 0.0;
    }

    const scalar k_B = plasma_->k_B_.value();
    const scalar E_ion = plasma_->ionizationEnergy_.value();
    const scalar m_e = plasma_->m_e_.value();

    const scalar h = 6.62607015e-34;
    const scalar pi = mathematical::pi;

    const scalar lambda_dB_sq = sqr(h) / (2.0 * pi * m_e * k_B * T);
    const scalar g_ratio = 2.0;

    const scalar K_saha = g_ratio * pow(2.0 / (lambda_dB_sq * n_atoms), 1.5)
                        * exp(-E_ion / (k_B * T));

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


Foam::tmp<Foam::volScalarField>
Foam::enhancedLIFTPhysics::calculateCharacteristicLength
(
    const volScalarField& alpha1
) const
{
    tmp<volScalarField> tL
    (
        new volScalarField
        (
            IOobject
            (
                "L_char",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            breakup_->minJetDiameter_
        )
    );

    volScalarField& L = tL.ref();

    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha) + dimensionedScalar("small", dimless/dimLength, SMALL));

    L = 1.0 / magGradAlpha;

    forAll(L, cellI)
    {
        L[cellI] = min
        (
            max(L[cellI], breakup_->minJetDiameter_.value()),
            breakup_->maxJetDiameter_.value()
        );
    }

    return tL;
}


void Foam::enhancedLIFTPhysics::calculateWeberNumber
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& sigma,
    const volScalarField& L_char,
    const volScalarField& alpha1
)
{
    const volScalarField magU(mag(U) + dimensionedScalar("smallU", dimVelocity, SMALL));
    const volScalarField sigmaSafe(sigma + dimensionedScalar("smallSigma", dimForce/dimLength, SMALL));

    // Calculate Weber number everywhere first
    breakup_->WeberNumber_ = rho * sqr(magU) * L_char / sigmaSafe;

    // Mask to zero outside interface region (0.01 < alpha < 0.99)
    // This ensures maximum Weber number is at the interface, not at boundaries
    scalarField& WeField = breakup_->WeberNumber_.primitiveFieldRef();
    const scalarField& alpha1Field = alpha1.primitiveField();

    forAll(WeField, cellI)
    {
        const scalar alpha = alpha1Field[cellI];
        if (alpha <= 0.01 || alpha >= 0.99)
        {
            WeField[cellI] = 0.0;
        }
    }

    breakup_->WeberNumber_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::enhancedLIFTPhysics::explosionFluxMultiplier
(
    const scalar T
) const
{
    if (!phaseExplosion_.valid() || !phaseExplosion_->enabled_ || T < phaseExplosion_->T_spinodal_.value())
    {
        return 1.0;
    }

    const scalar T_spin = phaseExplosion_->T_spinodal_.value();
    const scalar T_crit = phaseExplosion_->T_critical_.value();

    const scalar xi = min
    (
        (T - T_spin) / max(T_crit - T_spin, SMALL),
        1.0
    );

    return 1.0 + (phaseExplosion_->explosionMultiplier_ - 1.0) * sqr(xi);
}


void Foam::enhancedLIFTPhysics::updatePhaseExplosion
(
    const volScalarField& T,
    const volScalarField& alpha1,
    const volScalarField& rho
)
{
    if (!phaseExplosion_.valid() || !phaseExplosion_->enabled_)
    {
        return;
    }

    phaseExplosion_->explosiveMassSource_.primitiveFieldRef() = 0.0;
    phaseExplosion_->explosivePressure_.primitiveFieldRef() = 0.0;
    phaseExplosion_->explosionIndicator_.primitiveFieldRef() = 0.0;

    const scalar T_spin = phaseExplosion_->T_spinodal_.value();
    const scalar T_crit = phaseExplosion_->T_critical_.value();
    const scalar tau = phaseExplosion_->tau_explosion_.value();

    const scalarField& TField = T.primitiveField();
    const scalarField& alpha1Field = alpha1.primitiveField();
    const scalarField& rhoField = rho.primitiveField();

    scalarField& massSourceField = phaseExplosion_->explosiveMassSource_.primitiveFieldRef();
    scalarField& pressureField = phaseExplosion_->explosivePressure_.primitiveFieldRef();
    scalarField& indicatorField = phaseExplosion_->explosionIndicator_.primitiveFieldRef();

    label nExplosiveCells = 0;

    forAll(TField, cellI)
    {
        const scalar T_local = TField[cellI];
        const scalar alpha = alpha1Field[cellI];
        const scalar rho_local = rhoField[cellI];

        if (T_local > T_spin && alpha > 0.5)
        {
            const scalar xi = min
            (
                (T_local - T_spin) / max(T_crit - T_spin, SMALL),
                1.0
            );

            indicatorField[cellI] = xi;

            const scalar explosionRate = (rho_local / tau) * exp(10.0 * xi);
            massSourceField[cellI] = explosionRate;

            pressureField[cellI] = phaseExplosion_->p_critical_.value() * sqr(xi);

            nExplosiveCells++;
        }
    }

    phaseExplosion_->explosiveMassSource_.correctBoundaryConditions();
    phaseExplosion_->explosivePressure_.correctBoundaryConditions();
    phaseExplosion_->explosionIndicator_.correctBoundaryConditions();

    if (nExplosiveCells > 0 && Pstream::master() && mesh_.time().timeIndex() % 100 == 0)
    {
        Info<< "Phase explosion active:" << nl
            << "    Explosive cells = " << returnReduce(nExplosiveCells, sumOp<label>()) << nl
            << "    Max T = " << gMax(T) << " K" << nl
            << "    Max explosion pressure = " << gMax(phaseExplosion_->explosivePressure_)/1e6 << " MPa" << nl
            << "    Max mass source = " << gMax(phaseExplosion_->explosiveMassSource_) << " kg/m^3/s" << endl;
    }
}


void Foam::enhancedLIFTPhysics::updatePlasma
(
    const volScalarField& T,
    const volScalarField& rho,
    const volScalarField& alpha1
)
{
    if (!plasma_.valid() || !plasma_->enabled_)
    {
        return;
    }

    plasma_->ionizationDegree_.primitiveFieldRef() = 0.0;
    plasma_->n_electron_.primitiveFieldRef() = 0.0;
    plasma_->plasmaPressure_.primitiveFieldRef() = 0.0;
    plasma_->plasmaShielding_.primitiveFieldRef() = 0.0;

    const scalar T_ion = plasma_->T_ionization_.value();
    const scalar k_B = plasma_->k_B_.value();
    const scalar m_a = plasma_->m_atom_.value();

    const scalarField& TField = T.primitiveField();
    const scalarField& rhoField = rho.primitiveField();
    const scalarField& alpha1Field = alpha1.primitiveField();

    scalarField& alphaField = plasma_->ionizationDegree_.primitiveFieldRef();
    scalarField& n_e_Field = plasma_->n_electron_.primitiveFieldRef();
    scalarField& p_plasma_Field = plasma_->plasmaPressure_.primitiveFieldRef();
    scalarField& shieldField = plasma_->plasmaShielding_.primitiveFieldRef();

    label nPlasmaCells = 0;

    forAll(TField, cellI)
    {
        const scalar T_local = TField[cellI];
        const scalar rho_local = rhoField[cellI];
        const scalar alpha = alpha1Field[cellI];

        if (T_local > T_ion && alpha > 0.1)
        {
            const scalar n_local = max(rho_local / m_a, SMALL);
            const scalar alpha_ion = sahaIonization(T_local, n_local);

            alphaField[cellI] = alpha_ion;

            const scalar n_e = alpha_ion * n_local;
            n_e_Field[cellI] = n_e;

            const scalar p_plasma = 2.0 * n_e * k_B * T_local;
            p_plasma_Field[cellI] = p_plasma;

            const scalar n_crit = 1.1e27;
            const scalar shield = 1.0 - exp(-n_e / n_crit);
            shieldField[cellI] = shield;

            nPlasmaCells++;
        }
    }

    plasma_->ionizationDegree_.correctBoundaryConditions();
    plasma_->n_electron_.correctBoundaryConditions();
    plasma_->plasmaPressure_.correctBoundaryConditions();
    plasma_->plasmaShielding_.correctBoundaryConditions();

    if (nPlasmaCells > 0 && Pstream::master() && mesh_.time().timeIndex() % 100 == 0)
    {
        Info<< "Plasma formation active:" << nl
            << "    Plasma cells = " << returnReduce(nPlasmaCells, sumOp<label>()) << nl
            << "    Max ionization degree = " << gMax(plasma_->ionizationDegree_) << nl
            << "    Max electron density = " << gMax(plasma_->n_electron_) << " m^-3" << nl
            << "    Max plasma pressure = " << gMax(plasma_->plasmaPressure_)/1e6 << " MPa" << nl
            << "    Max laser shielding = " << gMax(plasma_->plasmaShielding_)*100.0 << " %" << endl;
    }
}


void Foam::enhancedLIFTPhysics::updateBreakup
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const volScalarField& rho
)
{
    if (!breakup_.valid() || !breakup_->enabled_)
    {
        return;
    }

    tmp<volScalarField> tsigma = mixture_.sigma();
    const volScalarField& sigma = tsigma();

    tmp<volScalarField> tL_char = calculateCharacteristicLength(alpha1);
    const volScalarField& L_char = tL_char();

    calculateWeberNumber(U, rho, sigma, L_char, alpha1);

    breakup_->breakupIndicator_.primitiveFieldRef() = 0.0;
    breakup_->dropletDiameter_.primitiveFieldRef() = 0.0;
    breakup_->breakupRate_.primitiveFieldRef() = 0.0;

    const scalarField& WeField = breakup_->WeberNumber_.primitiveField();
    const scalarField& alpha1Field = alpha1.primitiveField();
    const scalarField& LField = L_char.primitiveField();
    // Store tmp to avoid dangling reference
    tmp<scalarField> tUmagField = mag(U.primitiveField());
    const scalarField& UmagField = tUmagField();
    scalarField& indicatorField = breakup_->breakupIndicator_.primitiveFieldRef();
    scalarField& diameterField = breakup_->dropletDiameter_.primitiveFieldRef();
    scalarField& rateField = breakup_->breakupRate_.primitiveFieldRef();

    label nBreakupCells = 0;

    forAll(WeField, cellI)
    {
        const scalar We = WeField[cellI];
        const scalar alpha = alpha1Field[cellI];
        const scalar L = LField[cellI];
        const scalar Umag = UmagField[cellI];

        if (We > breakup_->We_critical_ && alpha > 0.01 && alpha < 0.99)
        {
            const scalar indicator = min((We - breakup_->We_critical_) / breakup_->We_critical_, 1.0);
            indicatorField[cellI] = indicator;

            const scalar d_drop = L / pow(We, 0.5);
            diameterField[cellI] = max(d_drop, breakup_->minJetDiameter_.value());

            const scalar tau_breakup = breakup_->breakupTimeCoeff_ * L / max(Umag, SMALL);
            rateField[cellI] = 1.0 / max(tau_breakup, SMALL);

            nBreakupCells++;
        }
    }

    breakup_->breakupIndicator_.correctBoundaryConditions();
    breakup_->dropletDiameter_.correctBoundaryConditions();
    breakup_->breakupRate_.correctBoundaryConditions();

    if (nBreakupCells > 0 && Pstream::master() && mesh_.time().timeIndex() % 100 == 0)
    {
        Info<< "Droplet breakup active:" << nl
            << "    Breakup cells = " << returnReduce(nBreakupCells, sumOp<label>()) << nl
            << "    Max Weber number = " << gMax(breakup_->WeberNumber_) << nl
            << "    Max breakup rate = " << gMax(breakup_->breakupRate_) << " s^-1" << nl
            << "    Min droplet diameter = " << gMin(breakup_->dropletDiameter_)*1e6 << " μm" << endl;
    }
}


void Foam::enhancedLIFTPhysics::applyBreakup
(
    volScalarField& alpha1,
    const volVectorField& U
)
{
    if (!breakup_.valid() || !breakup_->enabled_)
    {
        return;
    }

    const scalarField& indicatorField = breakup_->breakupIndicator_.primitiveField();
    scalarField& alpha1Field = alpha1.primitiveFieldRef();

    forAll(indicatorField, cellI)
    {
        const scalar indicator = indicatorField[cellI];

        if (indicator > 0.1)
        {
            const scalar smoothing = 0.01 * indicator;

            if (alpha1Field[cellI] > 0.5)
            {
                alpha1Field[cellI] = max(0.5, alpha1Field[cellI] - smoothing);
            }
            else if (alpha1Field[cellI] < 0.5)
            {
                alpha1Field[cellI] = min(0.5, alpha1Field[cellI] + smoothing);
            }
        }
    }

    alpha1.correctBoundaryConditions();
}


void Foam::enhancedLIFTPhysics::updateAll
(
    const volScalarField& T,
    const volScalarField& alpha1,
    const volScalarField& rho,
    const volVectorField& U
)
{
    updatePhaseExplosion(T, alpha1, rho);
    updatePlasma(T, rho, alpha1);
    updateBreakup(alpha1, U, rho);
}


// ************************************************************************* //
