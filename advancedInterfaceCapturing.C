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
    havePreviousRecoil_(false)
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

    pressureScale_ = dimensionedScalar
    (
        "pressureScale",
        dimPressure*dimTime,
        aicDict.lookupOrDefault<scalar>
        (
            "pressureScale",
            pressureScale_.value()
        )
    );
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
    // Check field validity using size and dimensions
    if (T_.size() != mesh_.nCells() || alpha1_.size() != mesh_.nCells())
    {
        FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
            << "Field size mismatch detected: T:" << T_.size() 
            << " alpha1:" << alpha1_.size()
            << " mesh:" << mesh_.nCells()
            << abort(FatalError);
    }

    // Ensure cached storage matches current mesh size
    if (previousRecoilPressure_.size() != mesh_.nCells())
    {
        previousRecoilPressure_.setSize(mesh_.nCells());
        previousRecoilPressure_ = 0.0;
        havePreviousRecoil_ = false;
    }

    // OPTIMIZED: Calculate once, reuse multiple times
    const scalar currentTime = mesh_.time().value();
    const bool master = Pstream::master();
    const scalarField& gasTField = T_.primitiveField();
    const scalarField& alpha1Field = alpha1_.primitiveField();
    const volScalarField& massRate = mixture_.dgdt();
    const scalarField* TlFieldPtr = nullptr;
    if (TlPtr_)
    {
        TlFieldPtr = &TlPtr_->primitiveField();
        if (TlFieldPtr->size() != gasTField.size())
        {
            FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
                << "Lattice temperature field size (" << TlFieldPtr->size()
                << ") does not match gas temperature size ("
                << gasTField.size() << ")"
                << abort(FatalError);
        }
    }

    if
    (
        massRate.size() != mesh_.nCells()
     || massRate.size() != gasTField.size()
     || massRate.size() != alpha1Field.size()
    )
    {
        FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
            << "Field size mismatch detected. massRate: " << massRate.size()
            << " T: " << gasTField.size()
            << " alpha1: " << alpha1Field.size()
            << " mesh: " << mesh_.nCells()
            << abort(FatalError);
    }

    const scalarField& massRateField = massRate.primitiveField();

    const scalar alphaWindow = Foam::max(alphaMax_ - alphaMin_, Foam::SMALL);
    const scalar recoilOnTemp = vaporTemp_.value() + recoilTempOffset_.value();
    const scalar evaporationOnTemp = vaporTemp_.value() + phaseChangeTempOffset_.value();    
    bool haveMetalCell = false;
    scalar maxTemp = 0.0;

    forAll(gasTField, cellI)
    {
        const scalar alpha = alpha1Field[cellI];
        const scalar gasT = gasTField[cellI];

        if (!std::isfinite(gasT))
        {
            FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
                << "Non-finite gas temperature value detected at cell " << cellI
                << ". Value: " << gasT
                << abort(FatalError);
        }

        if (TlFieldPtr)
        {
            const scalar latticeT = (*TlFieldPtr)[cellI];
            if (!std::isfinite(latticeT))
            {
                FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
                    << "Non-finite lattice temperature value detected at cell "
                    << cellI << ". Value: " << latticeT
                    << abort(FatalError);
            }
        }


        if (!std::isfinite(alpha))
        {
            FatalErrorIn("advancedInterfaceCapturing::calculateRecoilPressure()")
                << "Non-finite alpha1 value detected at cell " << cellI
                << ". Value: " << alpha
                << abort(FatalError);
        }
        scalar Tval = gasT;
        if (TlFieldPtr && alpha > 0.01)
        {
            Tval = (*TlFieldPtr)[cellI];
        }

        if (alpha > 0.01)
        {
            if (!haveMetalCell || Tval > maxTemp)
            {
                maxTemp = Tval;
            }
            haveMetalCell = true;
        }
    }
    const bool verbose = mesh_.time().controlDict().lookupOrDefault<Switch>("verbose", false);
    if (verbose && master)
    {
        Info<< "Calculating recoil pressure at t = " << currentTime
            << "s, max T = " << maxTemp
            << "K, recoil activation T = " << recoilOnTemp << "K"
            << ", evaporation T = " << evaporationOnTemp << "K" << endl;
    }
    // OPTIMIZED: Early exit using cached values
    if (!haveMetalCell || maxTemp < recoilOnTemp)
    {
        // OPTIMIZED: Use cached zero value instead of creating new one
        recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
        if (verbose && master)
        {
            Info<< "Temperature too low for recoil pressure" << endl;
        }
        previousRecoilPressure_ = recoilPressure_.primitiveField();
        havePreviousRecoil_ = true;
        recoilPressure_.correctBoundaryConditions();
        return;
    }
    // Initialize recoil pressure field
    recoilPressure_ = dimensionedScalar("zero", dimPressure, 0.0);
    // Constants for recoil pressure model sourced from dictionaries
    const scalar pressureScale = pressureScale_.value();
    const bool clampRecoil = clampRecoil_;
    const bool scaleRecoilMax = scaleRecoilMax_;
    const scalar recoilRelax = recoilRelax_;
    const bool applyRelaxation = recoilRelax < (1.0 - Foam::SMALL);
    const scalar massRateEps = 1e-12;
    // Access fields only once for efficiency and validate sizes
    scalarField& recoilField = recoilPressure_.primitiveFieldRef();
    
    const bool logRecoilSuppression =
        mesh_.time().controlDict().lookupOrDefault<Switch>
        (
            "logRecoilSuppression",
            false
        );
    const scalar invAlphaWindow = 1.0/alphaWindow;
    label suppressedCondensationCells = 0;

    // Compute recoil pressure based on evaporation rate
    forAll(gasTField, cellI)
    {
        const scalar alpha = alpha1Field[cellI];
        scalar alphaMask = (alpha - alphaMin_)*invAlphaWindow;
        alphaMask = Foam::min(Foam::max(alphaMask, scalar(0)), scalar(1));

        if (alphaMask <= SMALL || alphaMask >= (1.0 - SMALL))
        {
            if (massRateField[cellI] < -massRateEps)
            {
                ++suppressedCondensationCells;
            }
            continue;
        }

        scalar localTemp = gasTField[cellI];
        if (TlFieldPtr && alpha > 0.01)
        {
            localTemp = (*TlFieldPtr)[cellI];
        }

        if (localTemp < recoilOnTemp)
        {
            continue;
        }

        const scalar rawRate = massRateField[cellI];
        const scalar localRecoilMax = scaleRecoilMax
            ? recoilMax_ * alphaMask
            : recoilMax_;

        scalar localRecoil = 0.0;

        if (rawRate > massRateEps)
        {
            localRecoil = pressureScale * rawRate;
            if (clampRecoil)
            {
                localRecoil = Foam::min(localRecoil, localRecoilMax);
            }
        }
        else if (rawRate < -massRateEps)
        {
            ++suppressedCondensationCells;
            continue;
        }
        else
        {
            const scalar tempWindow = Foam::max(evaporationOnTemp - recoilOnTemp, Foam::SMALL);
            scalar tempWeight = (localTemp - recoilOnTemp)/tempWindow;
            tempWeight = Foam::min(Foam::max(tempWeight, scalar(0)), scalar(1));
            localRecoil = tempWeight * localRecoilMax;
        }

        localRecoil *= alphaMask;

        recoilField[cellI] = clampRecoil && !scaleRecoilMax
            ? Foam::min(localRecoil, recoilMax_)
            : localRecoil;
    }
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
    if (logRecoilSuppression && suppressedCondensationCells > 0 && master)
    {
        Info<< "advancedInterfaceCapturing: suppressed recoil in "
            << suppressedCondensationCells
            << " cells due to condensation mass flux." << endl;
    }
    // Ensure boundary conditions are correct
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
        Info<< "Phase-1 volume fraction = "
            << alpha1_.weightedAverage(mesh_.V()).value()
            << "  Min(alpha1) = " << min(alpha1_).value()
            << "  Max(alpha1) = " << max(alpha1_).value()
            << "  Max recoil pressure = " << max(recoilPressure_).value()
            << endl;
    }
}
void advancedInterfaceCapturing::write() const
{
    recoilPressure_.write();
}
} // End namespace Foam

