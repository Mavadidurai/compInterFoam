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

#include "dropletBreakupModel.H"
#include "twoPhaseMixtureThermo.H"
#include "fvc.H"
#include "Pstream.H"
#include "mathematicalConstants.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dropletBreakupModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dropletBreakupModel::dropletBreakupModel
(
    const fvMesh& mesh,
    const twoPhaseMixtureThermo& mixture,
    const dictionary& dict
)
:
    mesh_(mesh),
    mixture_(mixture),
    We_critical_
    (
        dict.lookupOrDefault<scalar>("We_critical", 10.0)
    ),
    minJetDiameter_
    (
        "minJetDiameter",
        dimLength,
        dict.lookupOrDefault<scalar>("minJetDiameter", 1e-7)  // 100 nm
    ),
    maxJetDiameter_
    (
        "maxJetDiameter",
        dimLength,
        dict.lookupOrDefault<scalar>("maxJetDiameter", 10e-6)  // 10 μm
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dropletBreakupModel::~dropletBreakupModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dropletBreakupModel::calculateCharacteristicLength
(
    const volScalarField& alpha1
) const
{
    // Estimate characteristic length from interface curvature
    // L ~ 1 / |∇α|

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
            minJetDiameter_
        )
    );

    volScalarField& L = tL.ref();

    // Calculate interface gradient magnitude
    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volScalarField magGradAlpha(mag(gradAlpha) + dimensionedScalar("small", dimless/dimLength, SMALL));

    // Characteristic length ~ 1 / |∇α|
    L = 1.0 / magGradAlpha;

    // Clamp to reasonable bounds
    forAll(L, cellI)
    {
        L[cellI] = min
        (
            max(L[cellI], minJetDiameter_.value()),
            maxJetDiameter_.value()
        );
    }

    return tL;
}


void Foam::dropletBreakupModel::calculateWeberNumber
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& sigma,
    const volScalarField& L_char
)
{
    // We = ρ * U^2 * L / σ

    const volScalarField magU(mag(U) + dimensionedScalar("smallU", dimVelocity, SMALL));
    const volScalarField sigmaSafe(sigma + dimensionedScalar("smallSigma", dimForce/dimLength, SMALL));

    WeberNumber_ = rho * sqr(magU) * L_char / sigmaSafe;

    WeberNumber_.correctBoundaryConditions();
}


void Foam::dropletBreakupModel::update
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const volScalarField& rho
)
{
    if (!enabled_)
    {
        return;
    }

    // Get surface tension
    tmp<volScalarField> tsigma = mixture_.sigma();
    const volScalarField& sigma = tsigma();

    // Calculate characteristic length
    tmp<volScalarField> tL_char = calculateCharacteristicLength(alpha1);
    const volScalarField& L_char = tL_char();

    // Calculate Weber number
    calculateWeberNumber(U, rho, sigma, L_char);

    // Reset breakup fields
    breakupIndicator_.primitiveFieldRef() = 0.0;
    dropletDiameter_.primitiveFieldRef() = 0.0;
    breakupRate_.primitiveFieldRef() = 0.0;

    const scalarField& WeField = WeberNumber_.primitiveField();
    const scalarField& alpha1Field = alpha1.primitiveField();
    const scalarField& LField = L_char.primitiveField();
    const scalarField& UmagField = mag(U.primitiveField())();

    scalarField& indicatorField = breakupIndicator_.primitiveFieldRef();
    scalarField& diameterField = dropletDiameter_.primitiveFieldRef();
    scalarField& rateField = breakupRate_.primitiveFieldRef();

    label nBreakupCells = 0;

    forAll(WeField, cellI)
    {
        const scalar We = WeField[cellI];
        const scalar alpha = alpha1Field[cellI];
        const scalar L = LField[cellI];
        const scalar Umag = UmagField[cellI];

        // Only consider interface cells (0.01 < alpha < 0.99) with high Weber number
        if (We > We_critical_ && alpha > 0.01 && alpha < 0.99)
        {
            // Breakup indicator: smooth transition above We_critical
            // indicator = (We - We_crit) / We_crit
            const scalar indicator = min((We - We_critical_) / We_critical_, 1.0);
            indicatorField[cellI] = indicator;

            // Estimated droplet diameter after breakup (Pilch & Erdman correlation)
            // d_drop ~ d_jet / We^0.5
            const scalar d_drop = L / pow(We, 0.5);
            diameterField[cellI] = max(d_drop, minJetDiameter_.value());

            // Breakup timescale: τ_breakup ~ d_jet / U
            const scalar tau_breakup = breakupTimeCoeff_ * L / max(Umag, SMALL);

            // Breakup rate: 1 / τ_breakup
            rateField[cellI] = 1.0 / max(tau_breakup, SMALL);

            nBreakupCells++;
        }
    }

    // Correct boundary conditions
    breakupIndicator_.correctBoundaryConditions();
    dropletDiameter_.correctBoundaryConditions();
    breakupRate_.correctBoundaryConditions();

    // Report statistics
    if (nBreakupCells > 0 && Pstream::master())
    {
        const scalar maxWe = gMax(WeberNumber_);
        const scalar maxBreakupRate = gMax(breakupRate_);
        const scalar minDropletD = gMin(dropletDiameter_);

        if (mesh_.time().timeIndex() % 100 == 0)
        {
            Info<< "Droplet breakup active:" << nl
                << "    Breakup cells = " << returnReduce(nBreakupCells, sumOp<label>()) << nl
                << "    Max Weber number = " << maxWe << nl
                << "    Max breakup rate = " << maxBreakupRate << " s^-1" << nl
                << "    Min droplet diameter = " << minDropletD*1e6 << " μm" << endl;
        }
    }
}


void Foam::dropletBreakupModel::applyBreakup
(
    volScalarField& alpha1,
    const volVectorField& U
)
{
    if (!enabled_)
    {
        return;
    }

    // Apply interface smoothing where breakup is active
    // This is a simplified model - a full implementation would use
    // interface remapping or Lagrangian particle tracking

    const scalarField& indicatorField = breakupIndicator_.primitiveField();
    scalarField& alpha1Field = alpha1.primitiveFieldRef();

    forAll(indicatorField, cellI)
    {
        const scalar indicator = indicatorField[cellI];

        if (indicator > 0.1)
        {
            // Smooth interface by slightly reducing sharp gradients
            // This represents sub-grid fragmentation
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


// ************************************************************************* //
