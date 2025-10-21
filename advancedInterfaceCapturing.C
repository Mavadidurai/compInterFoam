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
#include "IOstreams.H"
#include "token.H"
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
        dimensionedScalar
        (
            "pressureScale",
            dimPressure*dimTime,
            1.0
        )
    ),
    recoilMax_(0.0),
    recoilTempOffset_
    (
        dimensionedScalar("recoilTempOffset", dimTemperature, 0.0)
    ),
    phaseChangeTempOffset_
    (
        dimensionedScalar("phaseChangeTempOffset", dimTemperature, 0.0)
    ),
    clampRecoil_(false),
    scaleRecoilMax_(false),
    recoilRelax_(1.0),
    // Relaxed default bounds to reduce clipping of interface values
    alphaMin_(0.001),
    alphaMax_(0.999),
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
    momentumAccommodationCoeff_(0.18)
{
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

                if (value.dimensions() != defaultValue.dimensions())
                {
                    FatalIOErrorInFunction(aicDict)
                        << "Entry '" << key << "' has dimensions "
                        << value.dimensions()
                        << ", expected " << defaultValue.dimensions()
                        << exit(FatalIOError);
                }

                result = value;
            }

            return result;
        };

    meltingTemp_ = readTemperature("meltingTemperature", meltingTemp_);
    vaporTemp_ = readTemperature("vaporTemperature", vaporTemp_);
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
    const volScalarField* TlPtr = TlPtr_;
    if (!TlPtr)
    {
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        previousRecoilPressure_.setSize(mesh_.nCells());
        previousRecoilPressure_ = 0.0;
        havePreviousRecoil_ = false;
        recoilPressure_.correctBoundaryConditions();
        return;
    }

    const volScalarField& Tl = *TlPtr;
    const volScalarField& massFlux = mixture_.phaseChangeMassFlux();

    if
    (
        Tl.size() != mesh_.nCells()
     || massFlux.size() != mesh_.nCells()
     || alpha1_.size() != mesh_.nCells()
    )
    {
        FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
            << "Field size mismatch detected. Tl: " << Tl.size()
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

    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
    scalarField& recoilField = recoilPressure_.primitiveFieldRef();
    const scalarField& alpha1Field = alpha1_.primitiveField();
    const scalarField& TlField = Tl.primitiveField();
    const scalarField& massFluxField = massFlux.primitiveField();

    const scalar alphaWindow = Foam::max(alphaMax_ - alphaMin_, Foam::SMALL);
    const scalar invAlphaWindow = 1.0/alphaWindow;

    const scalar k_B = boltzmannConstant_.value();
    const scalar particleMass = vaporParticleMass_.value();
    const scalar betaMomentum = momentumAccommodationCoeff_;

    const bool verbose =
        mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);
    const bool master = Pstream::master();

    label clampCount = 0;
    scalar maxOvershootAmount = 0.0;
    scalar maxOvershootMag = 0.0;
    scalar maxOvershootLimit = 0.0;
    label maxOvershootCell = -1;

    forAll(recoilField, cellI)
    {
        const scalar alpha = alpha1Field[cellI];

        if (alpha < alphaMin_ || alpha > alphaMax_)
        {
            recoilField[cellI] = 0.0;
            continue;
        }

        const scalar jNet = massFluxField[cellI];
        if (Foam::mag(jNet) <= massRateEps_)
        {
            recoilField[cellI] = 0.0;
            continue;
        }

        const scalar Tlocal = Foam::max(TlField[cellI], scalar(0));
        const scalar vThermal = sqrt(k_B*Tlocal/particleMass);

        const scalar pRecoil = jNet*vThermal*(1.0 - betaMomentum);

        scalar alphaMask = (alpha - alphaMin_)*invAlphaWindow;
        alphaMask = Foam::min(Foam::max(alphaMask, scalar(0)), scalar(1));

        const scalar unclampedRecoil = pRecoil*alphaMask;

        // Honour the dictionary controlled clamp instead of a hard-coded ceiling.
        if (clampRecoil_)
        {
            scalar localMax = recoilMax_;
            if (scaleRecoilMax_)
            {
                localMax *= alphaMask;
            }

            localMax = Foam::max(localMax, scalar(0));

            const scalar unclampedMag = Foam::mag(unclampedRecoil);

            if (localMax > SMALL && unclampedMag > localMax)
            {
                ++clampCount;

                const scalar overshootAmount = unclampedMag - localMax;
                if (overshootAmount > maxOvershootAmount)
                {
                    maxOvershootAmount = overshootAmount;
                    maxOvershootMag = unclampedMag;
                    maxOvershootLimit = localMax;
                    maxOvershootCell = cellI;
                }

                recoilField[cellI] = Foam::min
                (
                    Foam::max(unclampedRecoil, -localMax),
                    localMax
                );
            }
            else
            {
                recoilField[cellI] = unclampedRecoil;
            }
        }
        else
        {
            recoilField[cellI] = unclampedRecoil;
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

        warn << ". Values were clamped." << endl;    }

    rampProgress_ = 1.0;

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

    const scalar recoilRelax = recoilRelax_;
    const bool applyRelaxation = recoilRelax < (1.0 - Foam::SMALL);
    if (applyRelaxation)
    {
        if (!havePreviousRecoil_)
        {
            previousRecoilPressure_ = recoilField;
            havePreviousRecoil_ = true;
        }
        else
        {
            forAll(recoilField, cellI)
            {
                recoilField[cellI] =
                    recoilRelax*recoilField[cellI]
                  + (1.0 - recoilRelax)*previousRecoilPressure_[cellI];
            }
            previousRecoilPressure_ = recoilField;
        }
    }
    else
    {
        previousRecoilPressure_ = recoilField;
        havePreviousRecoil_ = true;
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
            << alpha1Avg
            << "  Min(" << alpha1Name << ") = " << min(alpha1_).value()
            << "  Max(" << alpha1Name << ") = " << max(alpha1_).value()
            << "  Max recoil pressure = " << max(recoilPressure_).value()
            << endl;

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

