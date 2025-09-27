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

Class
    Foam::twoTemperatureModel

Description
    Two-temperature model for femtosecond laser-material interaction in LIFT
    processes. The model evolves separate electron and lattice temperatures,
    handles temperature-dependent properties, and tracks energy conservation.

SourceFiles
    twoTemperatureModel.C

    Dictionary entries:
        energyTolerance   Relative energy conservation tolerance used by
                          checkEnergyConservation
        minTe             Minimum temperature bound for Te and Tl [K]
        maxTe             Maximum temperature bound for Te and Tl [K]
        metalAmbientBlendWidth
                          Width of the metal/ambient blending layer used when
                          applying temperature bounds (optional, default 1e-3)
        klHighTThreshold  Lattice temperature [K] above which the high-
                          temperature conductivity correction activates
                          (optional, default 1000)
        klExponent        Exponent applied to (klHighTThreshold/Tl) when
                          reducing kl at high temperature (optional,
                          default 0.5, set to 0 to disable)
        gasMetalExchangeCoeff
                          Interfacial gas-metal coupling coefficient [W/m³/K].
                          Values in the 1e12-1e15 range provide moderate
                          coupling, while coefficients above 1e16 W/m³/K
                          enforce near-equilibrium interface temperatures and
                          can lead to stiff coupling.                     
\*---------------------------------------------------------------------------*/
#ifndef twoTemperatureModel_H
#define twoTemperatureModel_H
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "Function1.H"
namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class twoTemperatureModel Declaration
\*---------------------------------------------------------------------------*/
class twoTemperatureModel
{
private:
    // Private Data
        //- Reference to mesh
        const fvMesh& mesh_;
        //- Local copy of dictionary controls
        dictionary dict_;
        //- Metal phase indicator (e.g. alpha1)
        const volScalarField& metalFraction_;
        //- Ambient fallback temperature for non-metal cells
        dimensionedScalar ambientTemperature_;
        //- Temperature Fields
        volScalarField Te_;        //- Electron temperature [K]
        volScalarField Tl_;        //- Lattice temperature [K]
        //- Material Properties
        dimensionedScalar Ce_;     //- Electron heat capacity [J/m³/K]
        dimensionedScalar Cl_;     //- Lattice heat capacity [J/m³/K]
        dimensionedScalar G_;      //- Electron-phonon coupling [W/m³/K]
        dimensionedScalar De_;     //- Electron thermal diffusivity [m^2/s]
        dimensionedScalar gasMetalExchangeCoeff_;
        
        //- Optional temperature-dependent electron heat capacity Ce(T)
        autoPtr<Function1<scalar>> CeFunction_;
        //- Optional temperature-dependent coupling factor G(T)
        autoPtr<Function1<scalar>> GFunction_;
        //- Optional gas-metal exchange coefficient variation
        autoPtr<Function1<scalar>> gasMetalExchangeFunction_;        
        //- Energy Conservation Tracking
        mutable dimensionedScalar lastTotalEnergy_;  //- Previous total energy [J]
        mutable bool energyInitialized_;             //- Energy tracking flag
        //- Number of electron micro-steps executed in the previous solve
        label lastElectronSubCycles_;        
    // Private Member Functions
        //- Validate model parameters and dimensions
        bool validateParameters() const;
        //- Threshold for treating locations as metal-rich for validation
        scalar metalValidationThreshold() const;        
        //- Validate temperature field values
        bool validateFields() const;
        //- Active-metal mask used for bounds and coupling
        tmp<volScalarField> metalActiveMask(scalar cutoff) const;
        //- Metal indicator clamped for energy weighting
        tmp<volScalarField> clampedMetalFraction() const;        
        //- Apply temperature bounds using algebraic masking
        void applyTemperatureBounds
        (
            const volScalarField& activeMask,
            const dimensionedScalar& minTemp,
            const dimensionedScalar& maxTemp,
            const dimensionedScalar& ambient
        );
        //- Assemble, relax, solve, and correct the lattice equation
        void solveLatticeEquation
        (
            const volScalarField& metalEff,
            const volScalarField& G,
            const volScalarField& klField,
            const volScalarField& phaseChangeRelaxCoeff,
            const volScalarField& phaseChangeSource,
            const volScalarField& TlOld,
            const volScalarField& gasMetalHeatFlux,
            const dimensionedScalar& dtSub,
            const volScalarField& TlPrev
        );
        //- Assemble, relax, solve, and correct the electron equation
        void solveElectronEquation
        (
            const volScalarField& metalEff,
            const volScalarField& Ce,
            const volScalarField& ke,
            const volScalarField& G,
            const volScalarField& laserSource,
            const dimensionedScalar& dtSub,
            const volScalarField& TePrev
        );
        //- Integrated gas-metal coupling energy for this time-step
        dimensionedScalar couplingEnergy
        (
            const volScalarField& gasMetalHeatFluxMasked,
            const dimensionedScalar& dtDim
        ) const;
        //- Compute current total metal energy
        dimensionedScalar currentTotalEnergy() const;
        //- Emit optional energy diagnostics and tracking information
        void writeEnergyDiagnostics
        (
            const dimensionedScalar& laserEnergy,
            const dimensionedScalar& phaseChangeEnergy,
            const dimensionedScalar& couplingEnergy,
            const dimensionedScalar& electronEnergyBefore,
            const dimensionedScalar& latticeEnergyBefore,
            const Switch& energyDiagnostics
        ) const;
        //- Report temperature statistics when requested
        void writeSolveStatistics
        (
            const scalar residual,
            const Switch& statsSwitch
        ) const;
        //- Check energy conservation with expected energy change
        bool checkEnergyConservation
        (
            const dimensionedScalar& previousEnergy,
            const dimensionedScalar& expectedEnergyChange,
            const dimensionedScalar& currentEnergy,
            scalar& relativeError
        ) const;
        //- Detailed energy audit when requested
        void reportEnergyViolation
        (
            const dimensionedScalar& previousEnergy,
            const dimensionedScalar& totalEnergyInput,
            const dimensionedScalar& currentEnergy,
            scalar relativeError,
            scalar dt,
            const Switch& auditSwitch
        ) const;
        //- Update energy tracking variables
        void updateEnergyTracking() const;
        //- Update energy tracking with externally computed energy
        void updateEnergyTracking(const dimensionedScalar& currentEnergy) const;
        //- Determine the electron micro-step plan for a provided time-step
        label electronSubCycleCount(const dimensionedScalar& dt) const;
        //- Apply the micro-step plan and record diagnostics
        label activeElectronSubCycles(const dimensionedScalar& dt);        
public:
    // Constructors
        //- Construct from mesh and dictionary
        twoTemperatureModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const volScalarField& metalFraction
        );
        //- Destructor
        virtual ~twoTemperatureModel();
        //- Prevent copying
        twoTemperatureModel(const twoTemperatureModel&) = delete;
        void operator=(const twoTemperatureModel&) = delete;
    // Member Functions
        //- Solve temperature equations
        void solve
        (
            const volScalarField& laserSource,
            const volScalarField& phaseChangeSource,
            const volScalarField& phaseChangeRelaxCoeff,
            const volScalarField& gasMetalHeatFlux
        );
        //- Calculate electron thermal conductivity
        tmp<volScalarField> electronThermalConductivity() const;
        //- Calculate electron heat capacity
        tmp<volScalarField> electronHeatCapacity() const;
        //- Calculate electron-phonon coupling
        tmp<volScalarField> electronPhononCoupling() const;
        //- Gas-metal exchange coefficient for coupling with gas solver
        tmp<volScalarField> gasMetalExchangeCoeffField() const;        
        //- Check model validity
        bool valid() const;
        //- Write model information
        void write() const;
        //- Calculate thermal conductivities
        tmp<volScalarField> ke() const;  //- Electronic thermal conductivity
        tmp<volScalarField> kl() const;  //- Lattice thermal conductivity
        // Access Functions
        //- Return electron temperature (non-const)
        volScalarField& Te();        
        //- Return electron temperature
        const volScalarField& Te() const { return Te_; }
        //- Return lattice temperature (non-const)
        volScalarField& Tl();
        //- Return lattice temperature
        const volScalarField& Tl() const { return Tl_; }      
        //- Return electron heat capacity
        const dimensionedScalar& Ce() const { return Ce_; }
        //- Return lattice heat capacity
        const dimensionedScalar& Cl() const { return Cl_; }
        //- Return electron-phonon coupling
        const dimensionedScalar& G() const { return G_; }
        //- Return electron thermal diffusivity (uniform unless overridden)
        const dimensionedScalar& De() const { return De_; }
        //- Return the executed electron micro-steps from the previous solve
        label lastElectronSubCycles() const { return lastElectronSubCycles_; }
        //- Planned electron micro-steps for a supplied global time-step
        label plannedElectronSubCycles(const dimensionedScalar& dt) const;        
};
} // End namespace Foam
#endif
