/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Implementation of the advanced interface capturing helper for LIFT
    simulations. Provides recoil pressure calculation, optional clamping, and
    diagnostic output controlled via controlDict switches.
\*---------------------------------------------------------------------------*/
#include "advancedInterfaceCapturing.H"
#include "fvc.H"
#include "fvm.H"
#include "wallFvPatch.H"
#include "fvPatchField.H"
#include "Pstream.H"
#include "dimensionSet.H"
#include "mathematicalConstants.H"
#include "IOstreams.H"
#include "token.H"
#include "PstreamReduceOps.H"
#include <cmath>

namespace Foam
{
advancedInterfaceCapturing::advancedInterfaceCapturing
(
    const fvMesh& mesh,
    volScalarField& alpha1,
    const twoPhaseMixtureThermo& mixture,
    const volScalarField& T,
    const volScalarField* Tl
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    mixture_(mixture),
    T_(T),
    TlPtr_(Tl),
    meltingTemp_
    (
        dimensionedScalar
        (
            "meltingTemperature",
            dimTemperature,
            mixture.T_melt()
        )
    ),
    vaporTemp_
    (
        dimensionedScalar
        (
            "vaporTemperature",
            dimTemperature,
            mixture.T_vapor()
        )
    ),
    pressureScale_
    (
        // Stored with pressure dimensions for dictionary validation but used as
        // a dimensionless multiplier when computing recoil pressure.
        dimensionedScalar
        (
            "pressureScale",
            dimPressure,
            1.0
        )
    ),
    recoilMax_(1e8), // 100 MPa default cap, consistent with fs-LIFT measurements
    recoilTempOffset_
    (
        dimensionedScalar("recoilTempOffset", dimTemperature, 0.0)
    ),
    phaseChangeTempOffset_
    (
        dimensionedScalar("phaseChangeTempOffset", dimTemperature, 0.0)
    ),
    maxPhysicalTemperature_
    (
        dimensionedScalar("maxPhysicalTemperature", dimTemperature, GREAT)
    ),
    clampRecoil_(false),
    scaleRecoilMax_(false),
    recoilRelax_(1.0),
    // Relaxed default bounds to reduce clipping of interface values
    alphaMin_(0.001),
    alphaMax_(1.0),
    massRateEps_(1e-12),
    metalAlphaCutoff_(0.01),
    recoilPressure_
    (
        IOobject
        (
            "recoilPressure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
    previousRecoilPressure_(mesh.nCells(), 0.0),
    havePreviousRecoil_(false),
    rampProgress_(0.0),
    rampIncrement_(1.0),
    recoilMaxDelta_(0.0),
    recoilSmoothCoeff_(0.0),
    recoilSmoothIters_(0),
    femtosecondCase_
    (
        mesh.time().controlDict().lookupOrDefault<Switch>
        (
            "femtosecondCase",
            mesh.time().deltaTValue() <= 5e-12
        )
    ),
    boltzmannConstant_
    (
        dimensionedScalar
        (
            "boltzmannConstant",
            dimensionSet(1, 2, -2, -1, 0, 0, 0),
            1.38e-23
        )
    ),
    vaporParticleMass_
    (
        dimensionedScalar
        (
            "vaporParticleMass",
            dimensionSet(1, 0, 0, 0, 0, 0, 0),
            7.95e-26
        )
    ),
    momentumAccommodationCoeff_(1.0),   // Fix: was 0.18; 0.54*p_sat (Anisimov/Knight) already includes the Knudsen-layer back-flux, so 0.18 counted it twice
    recoilCoefficient_(0.54),
    relaxTimeIndex_(-1),
    pocketMass_(0),
    pocketDisp_(0),
    pocketRadius_(0),
    pocketPressure_(0),
    pocketTimeIndex_(-1),
    pocketTime_(mesh.time().value()),
    pocketTv_(0),
    pocketSavedMode_(),
    pocketSavedBinWidth_(0)
{
    // Confined vapour pocket: restore the 0-D state on a restart
    {
        IOdictionary pocketState
        (
            IOobject
            (
                "vapourPocket",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );
        pocketMass_ = pocketState.getOrDefault<scalar>("mass", 0);
        pocketDisp_ = pocketState.getOrDefault<scalar>("displacement", 0);
        pocketRadius_ = pocketState.getOrDefault<scalar>("radius", 0);
        pocketTv_ = pocketState.getOrDefault<scalar>("vapourTemperature", 0);
        pocketBinMass_ = pocketState.getOrDefault<scalarList>("binMass", scalarList());
        pocketBinDisp_ = pocketState.getOrDefault<scalarList>("binDisplacement", scalarList());
        pocketBinTv_ = pocketState.getOrDefault<scalarList>("binVapourTemperature", scalarList());
        pocketSavedMode_ = pocketState.getOrDefault<word>("mode", word::null);
        pocketSavedBinWidth_ = pocketState.getOrDefault<scalar>("binWidth", 0);
    }

    const dictionary& aicDict =
        mesh.time().controlDict().subOrEmptyDict("advancedInterfaceCapturing");

    const bool verbose =
        mesh.time().controlDict().lookupOrDefault<Switch>("verbose", false);

    auto readTemperature =
        [&](const word& key, const dimensionedScalar& defaultValue)
        {
            dimensionedScalar result(defaultValue);

            if (aicDict.found(key))
            {
                const dimensionedScalar value =
                    aicDict.lookupOrDefault<dimensionedScalar>(key, defaultValue);

                if (value.dimensions() == defaultValue.dimensions())
                {
                    result = value;
                }
                else if (value.dimensions() == dimless)
                {
                    result = dimensionedScalar
                    (
                        value.name(),
                        defaultValue.dimensions(),
                        value.value()
                    );
                }
                else
                {
                    FatalIOErrorInFunction(aicDict)
                        << "Entry '" << key << "' has dimensions "
                        << value.dimensions()
                        << ", expected " << defaultValue.dimensions()
                        << " or " << dimless
                        << exit(FatalIOError);
                }
            }

            return result;
        };

    meltingTemp_ = readTemperature("meltingTemperature", meltingTemp_);
    vaporTemp_ = readTemperature("vaporTemperature", vaporTemp_);
    const bool hasMaxPhysicalTemperature = aicDict.found("maxPhysicalTemperature");
    maxPhysicalTemperature_ = readTemperature
    (
        "maxPhysicalTemperature",
        maxPhysicalTemperature_
    );

    if
    (
        hasMaxPhysicalTemperature
     && maxPhysicalTemperature_.value() <= meltingTemp_.value()
    )
    {
        FatalErrorInFunction
            << "maxPhysicalTemperature (" << maxPhysicalTemperature_.value()
            << " K) must exceed the melting temperature ("
            << meltingTemp_.value() << " K). Values read from the"
            << " 'advancedInterfaceCapturing' dictionary."
            << abort(FatalError);
    }

    const bool master = Pstream::master();
    if (verbose && master)
    {
        Info<< "advancedInterfaceCapturing temperature bounds: melting = "
            << meltingTemp_.value() << " K, vapor = " << vaporTemp_.value()
            << " K" << endl;
    }

    if (meltingTemp_.value() >= vaporTemp_.value())
    {
        FatalErrorInFunction
            << "meltingTemperature (" << meltingTemp_.value()
            << " K) must be less than vaporTemperature ("
            << vaporTemp_.value() << " K). Values read from the"
            << " 'advancedInterfaceCapturing' dictionary."
            << abort(FatalError);
    }
    const scalar phaseChangeOffsetValue =
        aicDict.lookupOrDefault<scalar>
        (
            "phaseChangeTempOffset",
            phaseChangeTempOffset_.value()
        );

    phaseChangeTempOffset_ = dimensionedScalar
    (
        "phaseChangeTempOffset",
        dimTemperature,
        phaseChangeOffsetValue
    );

    const dimensionedScalar pressureDefault(pressureScale_);
    const dimensionedScalar rawPressureScale =
        aicDict.lookupOrDefault<dimensionedScalar>
        (
            "pressureScale",
            pressureDefault
        );

    dimensionedScalar pressureCfg(pressureDefault);

    if (rawPressureScale.dimensions() == pressureDefault.dimensions())
    {
        pressureCfg = rawPressureScale;
    }
    else if (rawPressureScale.dimensions() == dimless)
    {
        pressureCfg = dimensionedScalar
        (
            rawPressureScale.name(),
            pressureDefault.dimensions(),
            rawPressureScale.value()
        );
    }
    else
    {
        FatalIOErrorInFunction(aicDict)
            << "Entry 'pressureScale' has dimensions "
            << rawPressureScale.dimensions()
            << ", expected " << pressureDefault.dimensions()
            << " or " << dimless
            << exit(FatalIOError);
    }


    pressureScale_ = pressureCfg;
    recoilMax_ = aicDict.lookupOrDefault<scalar>("recoilMax", recoilMax_);

    // Honour the more common maxRecoilPressure control when present, and use
    // the tighter of the two caps.
    const scalar recoilMaxPressure =
        aicDict.lookupOrDefault<scalar>("maxRecoilPressure", recoilMax_);

    if (recoilMaxPressure > SMALL)
    {
        recoilMax_ =
            (recoilMax_ > SMALL)
          ? Foam::min(recoilMax_, recoilMaxPressure)
          : recoilMaxPressure;
    }
    
    const scalar recoilOffsetValue =
        aicDict.lookupOrDefault<scalar>
        (
            "recoilTempOffset",
            recoilTempOffset_.value()
        );

    recoilTempOffset_ = dimensionedScalar
    (
        "recoilTempOffset",
        dimTemperature,
        recoilOffsetValue
    );
    if (recoilTempOffset_.value() < 0)
    {
        FatalErrorInFunction
            << "recoilTempOffset (" << recoilTempOffset_.value()
            << ") must be non-negative"
            << abort(FatalError);
    }
    if (recoilTempOffset_.value() >= vaporTemp_.value())
    {
        FatalErrorInFunction
            << "recoilTempOffset (" << recoilTempOffset_.value()
            << ") must be less than vaporTemperature ("
            << vaporTemp_.value() << ")"
            << abort(FatalError);
    }
    if (recoilTempOffset_.value() > phaseChangeTempOffset_.value())
    {
        if (master)
        {
            WarningInFunction
                << "recoilTempOffset (" << recoilTempOffset_.value()
                << ") exceeds phaseChangeTempOffset ("
                << phaseChangeTempOffset_.value()
                << ") and will be limited." << endl;
        }
        recoilTempOffset_ = phaseChangeTempOffset_;
    }

    clampRecoil_ = aicDict.lookupOrDefault<Switch>
    (
        "clampRecoil",
        clampRecoil_
    );
    scaleRecoilMax_ = aicDict.lookupOrDefault<Switch>
    (
        "scaleRecoilMax",
        scaleRecoilMax_
    );
    if (clampRecoil_ && recoilMax_ <= SMALL)
    {
        FatalIOErrorInFunction(aicDict)
            << "recoilMax must be positive when clampRecoil is enabled"
            << exit(FatalIOError);
    }    
    if (aicDict.found("recoilRelax"))
    {
        recoilRelax_ = aicDict.lookupOrDefault<scalar>("recoilRelax", recoilRelax_);
    }
    else if (aicDict.found("relaxFactor"))
    {
        recoilRelax_ = aicDict.lookupOrDefault<scalar>("relaxFactor", recoilRelax_);
        if (master)
        {
            WarningInFunction
                << "Entry 'relaxFactor' in advancedInterfaceCapturing is deprecated. "
                << "Please use 'recoilRelax' instead." << endl;
        }
    }
    recoilRelax_ = Foam::min(Foam::max(recoilRelax_, scalar(0)), scalar(1));    
    alphaMin_ = aicDict.lookupOrDefault<scalar>("alphaMin", alphaMin_);
    alphaMax_ = aicDict.lookupOrDefault<scalar>("alphaMax", alphaMax_);
    massRateEps_ = Foam::max
    (
        aicDict.lookupOrDefault<scalar>("massRateEps", massRateEps_),
        scalar(0)
    );
    const scalar mixtureAlphaMin = mixture_.phaseChangeAlphaMin();
    if (alphaMin_ > mixtureAlphaMin)
    {
        if (master)
        {
            WarningInFunction
                << "alphaMin (" << alphaMin_
                << ") exceeds phase-change cutoff (" << mixtureAlphaMin
                << ") from twoPhaseMixtureThermo and will be limited." << endl;
        }

        alphaMin_ = mixtureAlphaMin;
    }
    if (femtosecondCase_)
    {
        // Keep the interface window wide open for femtosecond vaporization.
        alphaMin_ = 0.0;
        alphaMax_ = 1.0;
        phaseChangeTempOffset_ = dimensionedScalar
        (
            "phaseChangeTempOffset",
            dimTemperature,
            0.0
        );
        recoilTempOffset_ = dimensionedScalar
        (
            "recoilTempOffset",
            dimTemperature,
            0.0
        );
    }
    metalAlphaCutoff_ = aicDict.lookupOrDefault<scalar>
    (
        "metalAlphaCutoff",
        metalAlphaCutoff_
    );
    const label rampSteps = Foam::max
    (
        aicDict.lookupOrDefault<label>("recoilRampSteps", 1),
        label(1)
    );
    rampIncrement_ = 1.0/static_cast<scalar>(rampSteps);
    rampProgress_ = 0.0;
    recoilMaxDelta_ = Foam::max
    (
        aicDict.lookupOrDefault<scalar>("recoilMaxDelta", recoilMaxDelta_),
        scalar(0)
    );
    recoilSmoothCoeff_ = aicDict.lookupOrDefault<scalar>
    (
        "recoilSmoothCoeff",
        recoilSmoothCoeff_
    );
    recoilSmoothCoeff_ = Foam::min(Foam::max(recoilSmoothCoeff_, scalar(0)), scalar(1));
    recoilSmoothIters_ = Foam::max
    (
        aicDict.lookupOrDefault<label>("recoilSmoothIters", recoilSmoothIters_),
        label(0)
    );
    if (femtosecondCase_)
    {
        // Avoid damping impulses that occur within a few femtoseconds of pulse start.
        rampProgress_ = 1.0;
        rampIncrement_ = 1.0;
        recoilMaxDelta_ = 0.0;
        recoilSmoothCoeff_ = 0.0;
        recoilSmoothIters_ = 0;
        recoilMax_ = (recoilMax_ > SMALL) ? recoilMax_ : scalar(GREAT);
    }
    const dimensionedScalar defaultBoltzmann(boltzmannConstant_);
    const dimensionedScalar rawBoltzmann =
        aicDict.lookupOrDefault<dimensionedScalar>
        (
            "boltzmannConstant",
            defaultBoltzmann
        );

    if (rawBoltzmann.dimensions() != defaultBoltzmann.dimensions())
    {
        FatalIOErrorInFunction(aicDict)
            << "Entry 'boltzmannConstant' has dimensions "
            << rawBoltzmann.dimensions()
            << ", expected " << defaultBoltzmann.dimensions()
            << exit(FatalIOError);
    }
    if (rawBoltzmann.value() <= 0)
    {
        FatalErrorInFunction
            << "boltzmannConstant (" << rawBoltzmann.value()
            << ") must be positive"
            << abort(FatalError);
    }
    boltzmannConstant_ = rawBoltzmann;

    const dimensionedScalar defaultParticleMass(vaporParticleMass_);
    const dimensionedScalar rawParticleMass =
        aicDict.lookupOrDefault<dimensionedScalar>
        (
            "vaporParticleMass",
            defaultParticleMass
        );

    if (rawParticleMass.dimensions() != defaultParticleMass.dimensions())
    {
        FatalIOErrorInFunction(aicDict)
            << "Entry 'vaporParticleMass' has dimensions "
            << rawParticleMass.dimensions()
            << ", expected " << defaultParticleMass.dimensions()
            << exit(FatalIOError);
    }
    if (rawParticleMass.value() <= 0)
    {
        FatalErrorInFunction
            << "vaporParticleMass (" << rawParticleMass.value()
            << ") must be positive"
            << abort(FatalError);
    }
    vaporParticleMass_ = rawParticleMass;

    momentumAccommodationCoeff_ = aicDict.lookupOrDefault<scalar>
    (
        "momentumAccommodationCoeff",
        momentumAccommodationCoeff_
    );
    if
    (
        momentumAccommodationCoeff_ < scalar(0)
     || momentumAccommodationCoeff_ > scalar(1)
    )
    {
        FatalErrorInFunction
            << "momentumAccommodationCoeff (" << momentumAccommodationCoeff_
            << ") must lie in [0, 1]"
            << abort(FatalError);
    }

    recoilCoefficient_ = aicDict.lookupOrDefault<scalar>
    (
        "recoilCoefficient",
        recoilCoefficient_
    );
    if (recoilCoefficient_ <= 0 || recoilCoefficient_ > 1)
    {
        FatalErrorInFunction
            << "recoilCoefficient (" << recoilCoefficient_
            << ") must lie in (0, 1]; the Anisimov/Knight Knudsen-layer"
            << " analysis gives ~0.54"
            << abort(FatalError);
    }
    if (alphaMin_ >= alphaMax_)
    {
        FatalErrorInFunction
            << "alphaMin (" << alphaMin_ << ") must be less than alphaMax ("
            << alphaMax_ << ")" << abort(FatalError);
    }
    // Simple initialization, no calculations in constructor to avoid MPI issues
    if (verbose && master)
    {
        Info<< "Advanced interface capturing initialized" << endl;
    }
}
void advancedInterfaceCapturing::calculateRecoilPressure()
{
    const volScalarField* temperaturePtr = TlPtr_ ? TlPtr_ : &T_;
    if (!temperaturePtr)
    {
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        previousRecoilPressure_.setSize(mesh_.nCells());
        previousRecoilPressure_ = 0.0;
        havePreviousRecoil_ = false;
        rampProgress_ = 0.0;
        recoilPressure_.correctBoundaryConditions();
        return;
    }

    const volScalarField& temperature = *temperaturePtr;
    const volScalarField& massFlux = mixture_.phaseChangeMassFlux();

    const Time& time = mesh_.time();
    const scalar currentTime = time.value();
    const scalar laserEnd =
        time.controlDict().lookupOrDefault<scalar>("laserEndTime", GREAT);
    // Keep the recoil active for 50 ps after the end of the laser window
    const scalar recoilGrace = 50e-12;  // 50 picoseconds

    if (currentTime > laserEnd + recoilGrace)
    {
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        previousRecoilPressure_.setSize(mesh_.nCells());
        previousRecoilPressure_ = 0.0;
        havePreviousRecoil_ = false;
        rampProgress_ = 0.0;
        recoilPressure_.correctBoundaryConditions();
        return;
    }
    
    if
    (
        temperature.size() != mesh_.nCells()
     || massFlux.size() != mesh_.nCells()
     || alpha1_.size() != mesh_.nCells()
    )
    {
        FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
            << "Field size mismatch detected. temperature: " << temperature.size()
            << " massFlux: " << massFlux.size()
            << " alpha1: " << alpha1_.size()
            << " mesh: " << mesh_.nCells()
            << abort(FatalError);
    }

    if (previousRecoilPressure_.size() != mesh_.nCells())
    {
        previousRecoilPressure_.setSize(mesh_.nCells());
        previousRecoilPressure_ = 0.0;
        havePreviousRecoil_ = false;
    }

    // Freeze the relaxation/rate-limit reference once per time step so that
    // the repeated recoil updates inside one PIMPLE step (outer iterations,
    // pressure correctors, alpha sub-cycles) all damp against the SAME
    // previous-step field.  Updating the reference on every call made the
    // effective damping depend on the number of calls per step.
    const label currentTimeIndex = mesh_.time().timeIndex();
    if (currentTimeIndex != relaxTimeIndex_)
    {
        previousRecoilPressure_ = recoilPressure_.primitiveField();
        havePreviousRecoil_ = true;
        relaxTimeIndex_ = currentTimeIndex;
    }

    if (!femtosecondCase_)
    {
        if (rampProgress_ < scalar(1))
        {
            // Progressively enable recoil according to the configured ramp steps.
            rampProgress_ = Foam::min(rampProgress_ + rampIncrement_, scalar(1));
        }
        else
        {
            rampProgress_ = scalar(1);
        }
    }
    else
    {
        rampProgress_ = scalar(1);
    }

    const scalar rampFactor = rampProgress_;
    const bool applyRamp = (!femtosecondCase_)
        && rampFactor < (scalar(1) - Foam::SMALL);

    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
    scalarField& recoilField = recoilPressure_.primitiveFieldRef();
    const scalarField& alpha1Field = alpha1_.primitiveField();
    const scalarField& temperatureField = temperature.primitiveField();
    const scalarField& massFluxField = massFlux.primitiveField();

    const scalar maxPhysicalT = maxPhysicalTemperature_.value();
    
    const scalar alphaWindow = Foam::max(alphaMax_ - alphaMin_, Foam::SMALL);
    const scalar invAlphaWindow = 1.0/alphaWindow;

    const scalar k_B = boltzmannConstant_.value();
    const scalar particleMass = vaporParticleMass_.value();

    const scalar evaporationCoeff = mixture_.evaporationCoeff();
    const scalar safeParticleMass = Foam::max(Foam::mag(particleMass), VSMALL);
    const scalar R_specific = k_B/safeParticleMass;

    const bool validEvaporationCoeff = Foam::mag(evaporationCoeff) > Foam::SMALL;

    // With the Hertz-Knudsen flux j_net = alpha_e*p_sat/sqrt(2*pi*R*T),
    //   p_recoil = recoilCoefficient*p_sat
    //            = (recoilCoefficient/alpha_e)*j_net*sqrt(2*pi*R*T).
    // recoilCoefficient ~ 0.54 from the Knudsen-layer analysis (Anisimov
    // 1968; Knight, AIAA J. 17, 1979); the previous (2-beta_m)/2 prefactor
    // gave ~0.9*p_sat, ~70% above the accepted value.
    const scalar knightCoeff = validEvaporationCoeff
        ? recoilCoefficient_
           /Foam::max(Foam::mag(evaporationCoeff), Foam::SMALL)
        : 0.0;
    const scalar scaledKnightCoeff = pressureScale_.value()*knightCoeff;
    
    const bool verbose =
        mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);
    const bool master = Pstream::master();

    label clampCount = 0;
    scalar maxOvershootAmount = 0.0;
    scalar maxOvershootMag = 0.0;
    scalar maxOvershootLimit = 0.0;
    label maxOvershootCell = -1;

    const bool enforceUpper = alphaMax_ < (scalar(1) - Foam::SMALL);

    label interfaceCells = 0;
    label fluxCells = 0;
    label localHotInterfaceCells = 0;
    scalar localMaxMassFlux = 0.0;
    scalar localMaxRawMassFlux = 0.0;
    scalar localMaxActiveTemp = -GREAT;
    scalar localMinActiveTemp = GREAT;
    scalar localMaxInterfaceTemp = -GREAT;
    scalar localMinInterfaceTemp = GREAT;
    const scalar vaporThreshold =
        vaporTemp_.value() + phaseChangeTempOffset_.value();
    const scalar recoilThreshold =
        vaporTemp_.value() + recoilTempOffset_.value();

    // Option: confined vapour at the donor/carrier side. In the real
    // experiment the film's laser-side surface is bonded to the glass, so the
    // vapour formed there cannot escape: it pushes on the film with the full
    // saturation pressure instead of the 0.54*p_sat Knudsen-layer recoil of a
    // free surface. The model represents the glass as gas, so without this
    // both film surfaces evaporate freely and their recoils cancel once the
    // film is isothermal. confinedVapourFactor (default 1 =
    // off) multiplies the recoil on surface cells whose metal side faces away
    // from confinedVapourDirection (the unit vector pointing towards the glass).
    const dictionary& aicDictCV =
        mesh_.time().controlDict().subOrEmptyDict("advancedInterfaceCapturing");
    const scalar confinedVapourFactor =
        aicDictCV.lookupOrDefault<scalar>("confinedVapourFactor", 1.0);
    vector confinedDir =
        aicDictCV.lookupOrDefault<vector>("confinedVapourDirection", vector(0, 1, 0));
    confinedDir /= Foam::max(Foam::mag(confinedDir), VSMALL);
    tmp<volVectorField> tGradAlphaCV;
    if (Foam::mag(confinedVapourFactor - 1.0) > SMALL)
    {
        tGradAlphaCV = fvc::grad(alpha1_, "grad(alphaConfinedVapour)");
    }

    forAll(recoilField, cellI)
    {
        const scalar alpha = alpha1Field[cellI];

        if (alpha < alphaMin_ || (enforceUpper && alpha > alphaMax_))
        {
            recoilField[cellI] = 0.0;
            continue;
        }

        ++interfaceCells;

        if (alpha < metalAlphaCutoff_)
        {
            // Below the metal-dominated threshold: too little condensed
            // metal present for the kinetic recoil closure below (which
            // assumes emission from a metal surface) to be meaningful.
            recoilField[cellI] = 0.0;
            continue;
        }

        const scalar Tlocal = Foam::max(temperatureField[cellI], scalar(0));
        const scalar Tactive = Foam::min(Tlocal, maxPhysicalT);

        localMaxInterfaceTemp = Foam::max(localMaxInterfaceTemp, Tactive);
        localMinInterfaceTemp = Foam::min(localMinInterfaceTemp, Tactive);

        if (Tactive >= recoilThreshold)
        {
            ++localHotInterfaceCells;
        }

        if (Tactive < recoilThreshold)
        {
            recoilField[cellI] = 0.0;
            continue;
        }

        const scalar jNet = massFluxField[cellI];
        localMaxRawMassFlux = Foam::max(localMaxRawMassFlux, Foam::mag(jNet));
        if (Foam::mag(jNet) <= massRateEps_)
        {
            recoilField[cellI] = 0.0;
            continue;
        }

        ++fluxCells;

        localMaxActiveTemp = Foam::max(localMaxActiveTemp, Tactive);
        localMinActiveTemp = Foam::min(localMinActiveTemp, Tactive);
        localMaxMassFlux = Foam::max(localMaxMassFlux, Foam::mag(jNet));
        scalar sqrtTerm = 0.0;

        if (validEvaporationCoeff)
        {
            const scalar hkArgument =
                Foam::max
                (
                    Foam::constant::mathematical::twoPi*R_specific*Tactive,
                    scalar(0)
                );
            sqrtTerm = Foam::sqrt(hkArgument);
        }

        scalar alphaMask = scalar(1);
        if (scaleRecoilMax_)
        {
            scalar ramp = (alpha - alphaMin_)*invAlphaWindow;
            alphaMask = Foam::min(Foam::max(ramp, scalar(0)), scalar(1));
        }

        // Anisimov/Knight recoil closure:
        //   p_recoil = (recoilCoefficient/alpha_e) * j_net * sqrt(2*pi*R*T)
        //            = recoilCoefficient * p_sat(T)
        // The evaporation coefficient already scales j_net, so dividing by
        // alpha_e restores the saturation-pressure level. momentumAccommodationCoeff_
        // then scales the fraction of that kinetic momentum actually
        // transferred to the surface as recoil (documented in README.md).
        scalar pRecoil =
            momentumAccommodationCoeff_*scaledKnightCoeff*jNet*sqrtTerm;
        if (tGradAlphaCV.valid() && ((tGradAlphaCV()[cellI] & confinedDir) < 0))
        {
            pRecoil *= confinedVapourFactor;   // glass-side (confined) surface
        }
        const scalar rawRecoil = pRecoil;
        scalar recoilValue = rawRecoil;
        if (clampRecoil_)
        {
            scalar localMax = recoilMax_;
            if (scaleRecoilMax_)
            {
                localMax *= alphaMask;
            }

            localMax = Foam::max(localMax, scalar(0));

            if (localMax > SMALL)
            {
                const scalar rawMag = Foam::mag(rawRecoil);
                const scalar overshootAmount = Foam::max(rawMag - localMax, scalar(0));

                if (overshootAmount > maxOvershootAmount)
                {
                    maxOvershootAmount = overshootAmount;
                    maxOvershootMag = rawMag;
                    maxOvershootLimit = localMax;
                    maxOvershootCell = cellI;
                }
            }
        }

        // Honour the dictionary controlled clamp instead of a hard-coded ceiling.
        if (clampRecoil_)
        {
            scalar localMax = recoilMax_;
            if (scaleRecoilMax_)
            {
                localMax *= alphaMask;
            }

            localMax = Foam::max(localMax, scalar(0));

            const scalar unclampedMag = Foam::mag(rawRecoil);

            if (localMax > SMALL && unclampedMag > localMax)
            {
                ++clampCount;


                recoilValue = Foam::min
                (
                    Foam::max(rawRecoil, -localMax),
                    localMax
                );
            }
        }

        recoilField[cellI] = recoilValue;
    }
    
    if (applyRamp)
    {
        forAll(recoilField, cellI)
        {
            recoilField[cellI] *= rampFactor;
        }
    }
    
    if (clampCount > 0 && verbose && master)
    {
        Ostream& warn = WarningInFunction;
        warn    << clampCount
                << " recoil pressure value(s) exceeded the configured limit."
                << " Peak request " << maxOvershootMag/1e6
                << " MPa vs limit " << maxOvershootLimit/1e6
                << " MPa (overshoot " << maxOvershootAmount/1e6 << " MPa)";

        if (maxOvershootCell >= 0)
        {
            warn << " (e.g. cell " << maxOvershootCell << ')';
        }

        warn << ". Values were clamped." << endl;
    }

    if (recoilSmoothCoeff_ > SMALL && recoilSmoothIters_ > 0)
    {
        const labelListList& cellCells = mesh_.cellCells();
        const scalar coeff = Foam::min(Foam::max(recoilSmoothCoeff_, scalar(0)), scalar(1));
        scalarField workingField(recoilField);
        scalarField updatedField(recoilField.size(), 0.0);

        for (label iter = 0; iter < recoilSmoothIters_; ++iter)
        {
            forAll(workingField, cellI)
            {
                const labelList& neighbours = cellCells[cellI];
                if (neighbours.size())
                {
                    scalar sum = 0.0;
                    label count = 0;
                    forAll(neighbours, nbrI)
                    {
                        const label nbrCellI = neighbours[nbrI];
                        if (nbrCellI >= 0 && nbrCellI < workingField.size())
                        {
                            sum += workingField[nbrCellI];
                            ++count;
                        }
                    }
                    if (count > 0)
                    {
                        const scalar neighbourAvg = sum/static_cast<scalar>(count);
                        updatedField[cellI] =
                            (1.0 - coeff)*workingField[cellI] + coeff*neighbourAvg;
                        continue;
                    }
                }
                updatedField[cellI] = workingField[cellI];
            }
            workingField = updatedField;
        }

      recoilField = workingField;
    }
    scalar maxDeltaPerStep = 0.0;
    bool limitRecoilChange = false;
    if (recoilMaxDelta_ > SMALL)
    {
        const scalar deltaT = mesh_.time().deltaTValue();
        if (deltaT > SMALL)
        {
            // Convert the user-supplied change rate [Pa/s] into a per-step limit.
            maxDeltaPerStep = recoilMaxDelta_*deltaT;
            limitRecoilChange = maxDeltaPerStep > SMALL;
        }
    }

    if (limitRecoilChange)
    {
        forAll(recoilField, cellI)
        {
            const scalar prevValue = previousRecoilPressure_[cellI];
            const scalar delta = recoilField[cellI] - prevValue;

            if (Foam::mag(delta) > maxDeltaPerStep)
            {
                // Restrict how quickly the recoil can deviate from its previous value.
                const scalar limitedDelta =
                    (delta >= 0) ? maxDeltaPerStep : -maxDeltaPerStep;
                recoilField[cellI] = prevValue + limitedDelta;
            }
        }
    }
    const scalar recoilRelax = recoilRelax_;
    const bool applyRelaxation = recoilRelax < (1.0 - Foam::SMALL);
    if (applyRelaxation)
    {
        // Blend against the per-step frozen reference; the reference is
        // advanced once per time step at the top of this function.
        forAll(recoilField, cellI)
        {
            recoilField[cellI] =
                recoilRelax*recoilField[cellI]
              + (1.0 - recoilRelax)*previousRecoilPressure_[cellI];
        }
    }
    
    label globalInterfaceCells = interfaceCells;
    Foam::reduce(globalInterfaceCells, Foam::sumOp<label>());

    label globalFluxCells = fluxCells;
    Foam::reduce(globalFluxCells, Foam::sumOp<label>());

    label globalHotInterfaceCells = localHotInterfaceCells;
    Foam::reduce(globalHotInterfaceCells, Foam::sumOp<label>());

    scalar maxMassFlux = localMaxMassFlux;
    Foam::reduce(maxMassFlux, Foam::maxOp<scalar>());

    scalar maxRawMassFlux = localMaxRawMassFlux;
    Foam::reduce(maxRawMassFlux, Foam::maxOp<scalar>());

    scalar maxActiveTemp = localMaxActiveTemp;
    Foam::reduce(maxActiveTemp, Foam::maxOp<scalar>());

    scalar minActiveTemp = localMinActiveTemp;
    Foam::reduce(minActiveTemp, Foam::minOp<scalar>());

    scalar maxInterfaceTemp = localMaxInterfaceTemp;
    Foam::reduce(maxInterfaceTemp, Foam::maxOp<scalar>());

    scalar minInterfaceTemp = localMinInterfaceTemp;
    Foam::reduce(minInterfaceTemp, Foam::minOp<scalar>());

    const scalar maxRecoilMagFinal = Foam::gMax
    (
        Foam::mag(recoilPressure_.primitiveField())
    );

    if (verbose && master)
    {
        Info<< "Recoil diagnostics: "
            << globalFluxCells << " of " << globalInterfaceCells
            << " interface cells supplied mass flux above "
            << massRateEps_ << " kg/m^2/s" << nl;

        if (globalFluxCells > 0)
        {
            Info<< "  Max |j_net| = " << maxMassFlux
                << " kg/m^2/s";

            Info<< ", active temperature range = [" << minActiveTemp
                << ", " << maxActiveTemp << "] K" << nl;
        }
        else
        {
            Info<< "  No cells exceeded the mass flux epsilon ("
                << massRateEps_ << " kg/m^2/s)." << nl;

            if (globalInterfaceCells > 0)
            {
                const label globalFluxBelowEps =
                    globalInterfaceCells - globalFluxCells;

                if (globalFluxBelowEps > 0)
                {
                    Info<< "  " << globalFluxBelowEps
                        << " interface cell(s) were filtered by massRateEps." << nl;
                }

                Info<< "  Interface temperature range = ["
                    << minInterfaceTemp << ", " << maxInterfaceTemp
                    << "] K" << nl;

                if (maxRawMassFlux > massRateEps_)
                {
                    Info<< "  Max |j_net| = " << maxRawMassFlux
                        << " kg/m^2/s (filtered by massRateEps)." << nl;
                }

                if (globalHotInterfaceCells > 0)
                {
                    Info<< "  " << globalHotInterfaceCells
                        << " interface cell(s) exceeded the recoil threshold "
                        << recoilThreshold
                        << " K (phase-change threshold " << vaporThreshold
                        << " K) but still fell below massRateEps." << nl
                        << "    Inspect phaseChangeMassFlux, evaporationCoeff,"
                        << " and activation windows." << nl;
                }
                else
                {
                    Info<< "  Interface temperatures remained below the recoil"
                        << " threshold " << recoilThreshold << " K." << nl;
                }
            }
        }
        
        if (maxOvershootAmount > 0)
        {
            Info<< "  Recoil clamp requests occurred; peak raw request "
                << maxOvershootMag/1e6
                << " MPa vs limit " << maxOvershootLimit/1e6 << " MPa";

            if (maxOvershootCell >= 0)
            {
                Info<< " (e.g. cell " << maxOvershootCell << ')';
            }

            Info<< nl;
        }

        Info<< "  Max |recoilPressure| = " << maxRecoilMagFinal/1e6
            << " MPa" << endl;
    }

    // ------------------------------------------------------------------
    // Confined vapour pocket (approach A).
    // The glass side of the film is bonded to the donor glass, so the vapour
    // formed there is trapped between glass and film. Instead of the free-
    // surface Knudsen-layer recoil (0.54 p_sat) the glass-side surface feels
    // the pocket pressure of a single, uniform 0-D vapour volume:
    //   m_v   += dt * sum(j * |grad alpha| V)          (Hertz-Knudsen inflow)
    //   delta += dt * <u . (-n_glass)>_metal, r < r_p   (film moves away)
    //   V_p    = A_p (h0 + delta) + m_v / rho_metal     (gap + thinning)
    //   p_v    = min( <p_sat(T_s)>, m_v R T_s / V_p )   (saturated -> ideal gas)
    // p_v is put on every glass-side interface cell inside the pocket radius
    // and reaches the film through the same face traction as the recoil
    // (localMax; momentum balance verified 1.00). The compressibility of
    // the vapour lives in this ODE; the field solver stays incompressible.
    // ------------------------------------------------------------------
    {
        const dictionary pocketDict =
            aicDictCV.subOrEmptyDict("confinedVapourPocket");
        if (pocketDict.getOrDefault<Switch>("enabled", false))
        {
            // a uniform-mode run must not resume from a local-mode state
            if
            (
                pocketSavedMode_ == "local"
             && pocketDict.getOrDefault<word>("mode", "uniform") != "local"
            )
            {
                FatalErrorInFunction
                    << "confinedVapourPocket restart state was written in "
                    << "local mode but the case now uses mode "
                    << pocketDict.getOrDefault<word>("mode", "uniform")
                    << ". Restore mode local or remove uniform/vapourPocket."
                    << exit(FatalError);
            }
            const scalar h0 = pocketDict.getOrDefault<scalar>("minThickness", 1e-9);
            // only metal above this height (along confinedDir) is the
            // donor film; the receiver is metal too and diluted the film
            // velocity average ~100x (delta 0.3 nm instead of 25 nm)
            const scalar filmMinHeight =
                pocketDict.getOrDefault<scalar>("filmMinHeight", -GREAT);
            const scalar rhoMetal =
                pocketDict.getOrDefault<scalar>("metalDensity", 4515.0);
            const volVectorField gradA
            (
                fvc::grad(alpha1_, "grad(alphaConfinedVapour)")
            );
            const vectorField& gA = gradA.primitiveField();
            const volVectorField& Uf = mesh_.lookupObject<volVectorField>("U");
            const vectorField& Ui = Uf.primitiveField();
            const vectorField& Cc = mesh_.C().primitiveField();
            const scalarField& Vc = mesh_.V().field();
            // unmodified p_sat(T) published by the phase-change model; the
            // stored mass flux may be scaled or capped
            const volScalarField* pSatPub =
                mesh_.findObject<volScalarField>("phaseChangeSaturationPressure");
            auto radial = [&](const vector& x)
            {
                return Foam::mag(x - (x & confinedDir)*confinedDir);
            };
            auto glassSide = [&](label i)
            {
                return alpha1Field[i] >= alphaMin_ && alpha1Field[i] <= alphaMax_
                    && (gA[i] & confinedDir) < 0;
            };

            const word pocketMode =
                pocketDict.getOrDefault<word>("mode", "uniform");
            if (pocketMode == "local")
            {
            // ----------------------------------------------------------
            // Local mode: one vapour column per radial bin.
            // A uniform 0-D pocket over r_p (set by any evaporation above
            // massRateEps) pushed the film as a flat plug of radius 3.7 um
            // that sheared at the cold edge. The vapour cannot
            // equalise over micrometres on sub-ns times, so each bin keeps
            // its own vapour (no radial vapour flow):
            //   m_b   += dt * sum_b(j |grad a| V)
            //   d_b   += dt * <u.(-n_glass)>_metal,b
            //   V_b    = A_b (h0 + d_b) + m_b/rho_metal
            //   p_b    = min( <p_sat(T_s)>_b, m_b R T_b / V_b )
            //   m_b   <- min(m_b, p_sat V_b/(R T_b))   (excess condenses)
            // Glass-side cells get max(p_b, free-surface recoil).
            // ----------------------------------------------------------
            const scalar dr = pocketDict.getOrDefault<scalar>("binWidth", 1e-7);
            const scalar rMaxB = pocketDict.getOrDefault<scalar>("maxRadius", 1e-5);
            const label nB = Foam::max(label(1), label(rMaxB/dr + 0.5));
            // A restart must use the binning it was written with: a changed
            // binWidth/maxRadius would reassign or truncate the saved vapour
            // columns, a changed mode would ignore them
            if
            (
                !pocketBinMass_.empty()
             && (
                    (!pocketSavedMode_.empty() && pocketSavedMode_ != "local")
                 || (pocketSavedBinWidth_ > 0
                  && Foam::mag(pocketSavedBinWidth_ - dr) > 1e-9*dr)
                 || pocketBinMass_.size() != nB
                )
            )
            {
                FatalErrorInFunction
                    << "confinedVapourPocket restart state (mode "
                    << pocketSavedMode_ << ", binWidth " << pocketSavedBinWidth_
                    << ", " << pocketBinMass_.size() << " bins) does not match "
                    << "the current settings (mode local, binWidth " << dr
                    << ", " << nB << " bins). Restore the original settings or "
                    << "remove uniform/vapourPocket to restart the pocket."
                    << exit(FatalError);
            }
            if (pocketBinMass_.size() != nB)
            {
                pocketBinMass_.resize(nB, scalar(0));
                pocketBinDisp_.resize(nB, scalar(0));
                pocketBinTv_.resize(nB, scalar(0));
            }
            // cells at or beyond maxRadius belong to no bin (-1): they are
            // neither accumulated nor given a pocket pressure
            auto bin = [&](const vector& x) -> label
            {
                const scalar r = radial(x);
                return r < rMaxB ? Foam::min(nB - 1, label(r/dr)) : label(-1);
            };
            scalarList mdotB(nB, Zero), areaB(nB, Zero), TAB(nB, Zero),
                TjAB(nB, Zero), pSatAB(nB, Zero), mUB(nB, Zero), mmB(nB, Zero);
            forAll(Cc, i)
            {
                const label b = bin(Cc[i]);
                if (b < 0) continue;
                if
                (
                    alpha1Field[i] > 0.5
                 && (Cc[i] & confinedDir) > filmMinHeight
                )
                {
                    mUB[b] += alpha1Field[i]*Vc[i]*(-(Ui[i] & confinedDir));
                    mmB[b] += alpha1Field[i]*Vc[i];
                }
                if (!glassSide(i)) continue;
                const scalar dA = Foam::mag(gA[i])*Vc[i];
                const scalar T = Foam::max(temperatureField[i], scalar(1));
                const scalar j = Foam::max(massFluxField[i], scalar(0));
                areaB[b] += dA;
                mdotB[b] += j*dA;
                TAB[b] += T*j*dA;
                TjAB[b] += j*dA;
                // weighted by the evaporation flux j*dA (as T_v):
                // area weighting let cool partial cells pull the bin's p_sat
                // far below the hot surface's value
                if (pSatPub)
                {
                    pSatAB[b] += (*pSatPub)[i]*j*dA;
                }
                else if (validEvaporationCoeff)
                {
                    pSatAB[b] +=
                        j*Foam::sqrt(constant::mathematical::twoPi*R_specific*T)
                       /Foam::max(Foam::mag(evaporationCoeff), SMALL)*j*dA;
                }
            }
            Pstream::listCombineReduce(mdotB, plusEqOp<scalar>());
            Pstream::listCombineReduce(areaB, plusEqOp<scalar>());
            Pstream::listCombineReduce(TAB, plusEqOp<scalar>());
            Pstream::listCombineReduce(TjAB, plusEqOp<scalar>());
            Pstream::listCombineReduce(pSatAB, plusEqOp<scalar>());
            Pstream::listCombineReduce(mUB, plusEqOp<scalar>());
            Pstream::listCombineReduce(mmB, plusEqOp<scalar>());

            const scalar tNow = mesh_.time().value();
            const bool newStep = (tNow > pocketTime_*(1 + 1e-12) + VSMALL);
            const scalar dtStep = newStep ? tNow - pocketTime_ : 0;
            if (newStep)
            {
                pocketTime_ = tNow;
                pocketTimeIndex_ = currentTimeIndex;
            }
            scalarList pB(nB, Zero);
            scalar mTot = 0, pMax = 0, rAct = 0;
            forAll(pB, b)
            {
                if (newStep)
                {
                    pocketBinMass_[b] += mdotB[b]*dtStep;
                    pocketBinDisp_[b] +=
                        (mmB[b] > VSMALL ? mUB[b]/mmB[b] : 0)*dtStep;
                }
                if (TjAB[b] > VSMALL) pocketBinTv_[b] = TAB[b]/TjAB[b];
                const scalar Tv = pocketBinTv_[b];
                if (areaB[b] <= VSMALL || pocketBinMass_[b] <= 0 || Tv <= 0)
                {
                    continue;
                }
                const scalar Vb =
                    areaB[b]*(h0 + Foam::max(pocketBinDisp_[b], scalar(0)))
                  + pocketBinMass_[b]/Foam::max(rhoMetal, SMALL);
                const scalar pIdeal = pocketBinMass_[b]*R_specific*Tv/Vb;
                // The saturation cap (and condensation of the excess) needs
                // a current p_sat estimate, i.e. evaporation in this bin now.
                // Without one (e.g. the phase-change window has ended) the
                // stored vapour keeps its ideal-gas pressure.
                if (TjAB[b] > VSMALL)
                {
                    const scalar pSat = pSatAB[b]/TjAB[b];
                    pB[b] = Foam::min(pSat, pIdeal);
                    if (newStep && pIdeal > pSat)
                    {
                        pocketBinMass_[b] = pSat*Vb/(R_specific*Tv);
                    }
                }
                else
                {
                    pB[b] = pIdeal;
                }
                if (clampRecoil_) pB[b] = Foam::min(pB[b], recoilMax_);
                mTot += pocketBinMass_[b];
                pMax = Foam::max(pMax, pB[b]);
            }
            forAll(pB, b)
            {
                if (pB[b] > 0.01*pMax) rAct = (b + 1)*dr;
            }
            pocketPressure_ = pB[0];
            pocketMass_ = mTot;
            pocketDisp_ = pocketBinDisp_[0];
            pocketRadius_ = rAct;

            label nPocket = 0;
            forAll(Cc, i)
            {
                const label b = bin(Cc[i]);
                if (b >= 0 && glassSide(i) && alpha1Field[i] >= metalAlphaCutoff_)
                {
                    const scalar p = pB[b];
                    if (p > recoilField[i])
                    {
                        recoilField[i] = p;
                        ++nPocket;
                    }
                }
            }
            reduce(nPocket, sumOp<label>());

            if (newStep && master && (verbose || currentTimeIndex % 20 == 0))
            {
                // pressure at a few radii for the log
                auto pAt = [&](scalar r)
                {
                    return pB[Foam::min(nB - 1, label(r/dr))];
                };
                Info<< "    vapourPocket(local): t = " << tNow
                    << " m_v = " << mTot << " kg, p(0) = " << pB[0]
                    << " Pa, p(1um) = " << pAt(1e-6)
                    << " Pa, p(2um) = " << pAt(2e-6)
                    << " Pa, p(3um) = " << pAt(3e-6)
                    << " Pa, p_max = " << pMax
                    << " Pa, r(p>1% max) = " << rAct
                    << " m, d(0) = " << pocketBinDisp_[0]
                    << " m, T_v(0) = " << pocketBinTv_[0]
                    << " K, cells raised = " << nPocket << endl;
            }
            if (mesh_.time().writeTime() && !mesh_.time().subCycling())
            {
                IOdictionary pocketState
                (
                    IOobject
                    (
                        "vapourPocket",
                        mesh_.time().timeName(),
                        "uniform",
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                );
                pocketState.add("mode", pocketMode);
                pocketState.add("binWidth", dr);
                pocketState.add("binMass", pocketBinMass_);
                pocketState.add("binDisplacement", pocketBinDisp_);
                pocketState.add("binVapourTemperature", pocketBinTv_);
                // ASCII: the bin lists are read back through a dictionary token
                // stream, which cannot read binary list blobs
                pocketState.regIOobject::writeObject
                (
                    IOstreamOption(IOstreamOption::ASCII),
                    true
                );
            }
            }
            else
            {
            // pass 1: evaporation on the glass side
            scalar mdot = 0, rEvap = 0;
            forAll(Cc, i)
            {
                if (!glassSide(i)) continue;
                const scalar j = massFluxField[i];
                if (j > massRateEps_)
                {
                    mdot += j*Foam::mag(gA[i])*Vc[i];
                    rEvap = Foam::max(rEvap, radial(Cc[i]));
                }
            }
            reduce(mdot, sumOp<scalar>());
            reduce(rEvap, maxOp<scalar>());
            const scalar rP = Foam::max(pocketRadius_, rEvap);

            // pass 2: pocket area, surface temperature, saturation cap, film velocity
            scalar area = 0, TA = 0, TjA = 0, pSatA = 0, evapA = 0, mU = 0, mm = 0;
            forAll(Cc, i)
            {
                const scalar r = radial(Cc[i]);
                if (r > rP) continue;
                if
                (
                    alpha1Field[i] > 0.5
                 && (Cc[i] & confinedDir) > filmMinHeight
                )
                {
                    mU += alpha1Field[i]*Vc[i]*(-(Ui[i] & confinedDir));
                    mm += alpha1Field[i]*Vc[i];
                }
                if (!glassSide(i)) continue;
                const scalar dA = Foam::mag(gA[i])*Vc[i];
                const scalar T = Foam::max(temperatureField[i], scalar(1));
                area += dA;
                const scalar j = massFluxField[i];
                if (j > massRateEps_ && (pSatPub || validEvaporationCoeff))
                {
                    // vapour temperature: evaporation-flux weighted (the
                    // vapour is produced mostly at the hot centre)
                    TA += T*j*dA;
                    TjA += j*dA;
                    const scalar pSat =
                        pSatPub
                      ? (*pSatPub)[i]
                      : j*Foam::sqrt(constant::mathematical::twoPi*R_specific*T)
                       /Foam::max(Foam::mag(evaporationCoeff), SMALL);
                    pSatA += pSat*dA;
                    evapA += dA;
                }
            }
            reduce(area, sumOp<scalar>());
            reduce(TA, sumOp<scalar>());
            reduce(TjA, sumOp<scalar>());
            reduce(pSatA, sumOp<scalar>());
            reduce(evapA, sumOp<scalar>());
            reduce(mU, sumOp<scalar>());
            reduce(mm, sumOp<scalar>());

            // advance the 0-D state by the simulation time elapsed since the
            // last update (robust to alpha sub-cycles and repeated calls)
            const scalar tNow = mesh_.time().value();
            const bool newStep = (tNow > pocketTime_*(1 + 1e-12) + VSMALL);
            if (newStep)
            {
                const scalar dtStep = tNow - pocketTime_;
                pocketMass_ += mdot*dtStep;
                pocketDisp_ += (mm > VSMALL ? mU/mm : 0)*dtStep;
                pocketRadius_ = rP;
                pocketTime_ = tNow;
                pocketTimeIndex_ = currentTimeIndex;
            }

            // vapour temperature: flux-weighted while evaporating; afterwards
            // keep the last value (the pocket is not cooled in this 0-D model)
            if (TjA > VSMALL) pocketTv_ = TA/TjA;
            const scalar Ts = pocketTv_;
            const scalar Vp =
                area*(h0 + Foam::max(pocketDisp_, scalar(0)))
              + pocketMass_/Foam::max(rhoMetal, SMALL);
            const scalar pIdeal =
                Vp > VSMALL ? pocketMass_*R_specific*Ts/Vp : 0;
            const scalar pSatMean = evapA > VSMALL ? pSatA/evapA : 0;
            pocketPressure_ =
                pSatMean > 0 ? Foam::min(pSatMean, pIdeal) : pIdeal;
            // honour the configured recoil cap for the pocket too
            if (clampRecoil_)
            {
                pocketPressure_ = Foam::min(pocketPressure_, recoilMax_);
            }

            label nPocket = 0;
            forAll(Cc, i)
            {
                if
                (
                    glassSide(i)
                 && alpha1Field[i] >= metalAlphaCutoff_
                 && radial(Cc[i]) <= rP
                )
                {
                    recoilField[i] = pocketPressure_;
                    ++nPocket;
                }
            }
            reduce(nPocket, sumOp<label>());

            if (newStep)
            {
                if (master && (verbose || currentTimeIndex % 20 == 0))
                {
                    Info<< "    vapourPocket: t = " << mesh_.time().value()
                        << " m_v = " << pocketMass_ << " kg, delta = "
                        << pocketDisp_ << " m, r_p = " << rP
                        << " m, A_p = " << area << " m2, V_p = " << Vp
                        << " m3, T_s = " << Ts << " K, p_sat = " << pSatMean
                        << " Pa, p_ideal = " << pIdeal << " Pa -> p_v = "
                        << pocketPressure_ << " Pa on " << nPocket
                        << " cells" << endl;
                }
            }
            // restart state: written on the real write step only (alpha
            // sub-cycles report intermediate times with writeTime() true),
            // on every call of that step so the final state is stored
            {
                if (mesh_.time().writeTime() && !mesh_.time().subCycling())
                {
                    IOdictionary pocketState
                    (
                        IOobject
                        (
                            "vapourPocket",
                            mesh_.time().timeName(),
                            "uniform",
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    );
                    pocketState.add("mass", pocketMass_);
                    pocketState.add("displacement", pocketDisp_);
                    pocketState.add("radius", pocketRadius_);
                    pocketState.add("pressure", pocketPressure_);
                    pocketState.add("vapourTemperature", pocketTv_);
                    pocketState.regIOobject::write();
                }
            }
            }   // uniform mode
        }
    }

    recoilPressure_.correctBoundaryConditions();
}
void Foam::advancedInterfaceCapturing::correct()
{
    const bool verbose = mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);
    const bool master = Pstream::master();
    if (verbose && master)
    {
        Info<< "Performing simplified interface capturing" << endl;
    }
    // Update the recoil pressure field every invocation so it follows fast transients.
    calculateRecoilPressure();

    alpha1_.correctBoundaryConditions();
    if (verbose && master)
    {
        const word& alpha1Name = alpha1_.name();
        const scalar alpha1Avg = alpha1_.weightedAverage(mesh_.V()).value();

        Info<< "Phase-1 volume fraction (" << alpha1Name << ") = "
            << alpha1Avg << nl
            << "  (domain average - includes donor and receiver substrates)"
            << nl
            << "  Min(" << alpha1Name << ") = " << min(alpha1_).value() << nl
            << "  Max(" << alpha1Name << ") = " << max(alpha1_).value() << nl
            << "  Max recoil pressure = " << max(recoilPressure_).value()/1e6
            << " MPa" << endl;

        // The domain-average alpha is dominated by the stationary substrates,
        // so it barely changes even when the thin donor film evolves. Report
        // an additional diagnostic restricted to cells that currently
        // participate in phase change so users can track the active film.
        const volScalarField& massFlux = mixture_.phaseChangeMassFlux();
        const scalarField& alphaField = alpha1_.internalField();
        const scalarField& volField = mesh_.V();
        const scalarField& massFluxField = massFlux.internalField();

        scalar activeVolume = 0.0;
        scalar activeMetalVolume = 0.0;

        forAll(massFluxField, cellI)
        {
            if (mag(massFluxField[cellI]) > massRateEps_)
            {
                const scalar cellVol = volField[cellI];
                activeVolume += cellVol;
                activeMetalVolume += alphaField[cellI]*cellVol;
            }
        }

        reduce(activeVolume, sumOp<scalar>());
        reduce(activeMetalVolume, sumOp<scalar>());

        if (activeVolume > VSMALL)
        {
            const scalar activeAlphaAvg = activeMetalVolume/activeVolume;

            Info<< "  Active-film alpha average (phase-change cells) = "
                << activeAlphaAvg
                << " (evaluated over " << activeVolume*1e18
                << " µm^3 participating in mass transfer)" << endl;
        }
        else
        {
            Info<< "  Active-film alpha average (phase-change cells) = n/a"
                << " (no cells above massRateEps)" << endl;
        }

        const volScalarField* alpha2Ptr = nullptr;

        if (mesh_.foundObject<volScalarField>("alpha.air"))
        {
            alpha2Ptr = &mesh_.lookupObject<volScalarField>("alpha.air");
        }
        else if (mesh_.foundObject<volScalarField>("alpha2"))
        {
            alpha2Ptr = &mesh_.lookupObject<volScalarField>("alpha2");
        }

        if (alpha2Ptr)
        {
            const volScalarField& alpha2Field = *alpha2Ptr;
            const word& alpha2Name = alpha2Field.name();

            Info<< "Phase-2 volume fraction (" << alpha2Name << ") = "
                << alpha2Field.weightedAverage(mesh_.V()).value()
                << "  Min(" << alpha2Name << ") = " << min(alpha2Field).value()
                << "  Max(" << alpha2Name << ") = " << max(alpha2Field).value()
                << endl;
        }
        else
        {
            const scalar alpha2Avg = scalar(1) - alpha1Avg;

            Info<< "Phase-2 volume fraction = "
                << alpha2Avg
                << "  Min(1-" << alpha1Name << ") = "
                << min(scalar(1) - alpha1_).value()
                << "  Max(1-" << alpha1Name << ") = "
                << max(scalar(1) - alpha1_).value()
                << endl;
        }
    }
}
void advancedInterfaceCapturing::write() const
{
    recoilPressure_.write();
}
} // End namespace Foam

