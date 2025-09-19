/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration    |
    \\  /    A nd          | www.openfoam.com
     \\/     M anipulation |
-------------------------------------------------------------------------------
    Description
    Two-temperature model for femtosecond laser-material interaction in LIFT process.
    Models separate electron and lattice temperatures, including:
    - Electron-phonon coupling
    - Temperature-dependent thermal properties
    - Energy conservation tracking
    - Laser energy absorption
    The model solves coupled equations for:
    - Electron temperature (Te)
    - Lattice temperature (Tl)
    With temperature-dependent properties:
    - Electron thermal conductivity (ke)
    - Lattice thermal conductivity (kl)
    - Electron-phonon coupling coefficient (G)   
SourceFiles
    twoTemperatureModel.C
        Dictionary entries:
        energyTolerance   Relative energy conservation tolerance
                          used by checkEnergyConservation
        minTe             Minimum temperature bound for Te and Tl [K]
        maxTe             Maximum temperature bound for Te and Tl [K]
        klHighTThreshold  Lattice temperature [K] above which the high-
                          temperature conductivity correction activates
                          (optional, default 1000)
        klExponent        Exponent applied to (klHighTThreshold/Tl) when
                          reducing kl at high temperature (optional,
                          default 0.5, set to 0 to disable)
        gasMetalExchangeCoeff  Interfacial gas-metal coupling coefficient
                          [W/m³/K]. Values in the 1e12-1e15 range provide
                          moderate coupling, while coefficients above 1e16
                          W/m³/K enforce near-equilibrium interface
                          temperatures and can lead to stiff coupling.                          
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
        //- Reference to dictionary
        const dictionary& dict_;
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
    // Private Member Functions
        //- Validate model parameters and dimensions
        bool validateParameters() const;
        //- Threshold for treating locations as metal-rich for validation
        scalar metalValidationThreshold() const;        
        //- Validate temperature field values
        bool validateFields() const;
        //- Abort if Te/Tl develop invalid values during the solve
        void guardTemperatureFields
        (
            const char* context,
            label sweep
        ) const;        
        //- Check energy conservation with expected energy change
        bool checkEnergyConservation
        (
            const dimensionedScalar& expectedEnergyChange
        ) const;
        //- Update energy tracking variables
        void updateEnergyTracking() const;
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
};
} // End namespace Foam
#endif
