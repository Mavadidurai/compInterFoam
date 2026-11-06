#ifndef twoPhaseMixtureThermo_H
#define twoPhaseMixtureThermo_H

#include "rhoThermo.H"
#include "psiThermo.H"
#include "twoPhaseMixture.H"
#include "interfaceProperties.H"
#include "dimensionedScalar.H"
#include "error.H"
#include "Tuple2.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class twoPhaseMixtureThermo Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMixtureThermo
:
    public psiThermo,
    public twoPhaseMixture,
    public interfaceProperties
{
    // Private data

        //- Thermo-package of phase 1
        autoPtr<rhoThermo> thermo1_;

        //- Thermo-package of phase 2
        autoPtr<rhoThermo> thermo2_;

        //- Phase-change/evaporation latent heat used for mass-transfer (J/kg)
        scalar latentHeat_;

        //- Solidus temperature (K)
        scalar T_melt_;

        //- Vaporisation temperature (K)
        scalar T_vapor_;
        
        //- Reference saturation pressure at T_vapor (Pa)
        scalar p_ref_;
        
        //- Optional physical temperature ceiling for phase-change physics (K)
        scalar maxPhaseChangeTemperature_;
        
        //- Minimum permissible time-step for phase change (s)
        mutable scalar dtFloor_;

        //- Ideal-gas constant used in Hertz-Knudsen relation [J/(kg K)]
        //  Must be supplied via controlDict.phaseChangeCoeffs
        scalar gasConstant_;

        //- Accommodation coefficient for evaporation/condensation []
        scalar evaporationCoeff_;

        //- Relaxation time applied to the explicit source term [s]
        scalar relaxationTime_;
        //- Maximum volumetric phase-change source [W/m^3]
        dimensionedScalar maxPhaseChangeSource_;

        //- Phase-fraction window for evaluating evaporation []
        scalar alphaMin_;
        scalar alphaMax_;

        //- Minimum |grad(alpha)|*cellSize identifying a resolved interface;
        //  cells below it carry no Hertz-Knudsen source (bulk metal) []
        scalar interfaceGradientCutoff_;

        //- Restrict phase change to temperatures above T_vapor
        Switch onlyAboveVapor_;

        //- Optional activation windows for phase change [s]
        List<Tuple2<scalar, scalar>> activationWindows_;
        //- Laser volumetric heat source (W/m^3)
        volScalarField Q_laser_;

        //- Phase change source term (explicit part)
        volScalarField phaseChangeSource_;

        //- Relaxation coefficient used for implicit phase-change treatment
        volScalarField phaseChangeRelaxCoeff_;
        //- Mass-transfer rate due to phase change [1/s]
        volScalarField dgdt_;
        //- Physical phase-change mass flux [kg/(m^2 s)]
        volScalarField phaseChangeMassFlux_;
        //- Volumetric mass-transfer rate [1/s] consistent with the
        //  unconditional Hertz-Knudsen phaseChangeSource_/phaseChangeMassFlux_
        //  computed in computePhaseChange(). Used as the default dgdt_ when
        //  the optional massTransferCoeffs dict is not configured, so latent
        //  heat removed from the lattice always has matching mass leaving
        //  alpha1 (see computeMassTransfer()).
        volScalarField phaseChangeVolumetricRate_;
        //- Kinematic viscosity of each phase from transportProperties
        dimensionedScalar nu1_;
        dimensionedScalar nu2_;

        //- Density of each phase from transportProperties
        dimensionedScalar rho1_;
        dimensionedScalar rho2_;
        //- Cached lattice heat capacity from the two-temperature model [J/m^3/K]
        //  Defaults to 2.5e6 J/m^3/K if twoTemperatureProperties are absent
        dimensionedScalar ClTTM_;

        //- Surface tension model
        mutable autoPtr<surfaceTensionModel> sigmaModel_;

        //- Refresh cached two-temperature properties from the control dictionary
        void updateTwoTemperatureCache();
public:

    //- Runtime type information
    TypeName("twoPhaseMixtureThermo");


    // Constructors

        //- Construct from components
        twoPhaseMixtureThermo
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~twoPhaseMixtureThermo();


    // Member Functions

        const rhoThermo& thermo1() const
        {
            return *thermo1_;
        }

        const rhoThermo& thermo2() const
        {
            return *thermo2_;
        }

        rhoThermo& thermo1()
        {
            return *thermo1_;
        }

        rhoThermo& thermo2()
       {
            return *thermo2_;
        }
       //- Return latent heat
        dimensionedScalar latentHeat() const;

        //- Return melting temperature
        scalar T_melt() const
        {
            return T_melt_;
        }
        //- Return vaporisation temperature
        scalar T_vapor() const
        {
            return T_vapor_;
        }
        
        //- Return the minimum alpha threshold for phase change
        scalar phaseChangeAlphaMin() const;
        //- Access kinematic viscosities read from transportProperties
        const dimensionedScalar& nu1() const { return nu1_; }
        const dimensionedScalar& nu2() const { return nu2_; }

        //- Access densities read from transportProperties
        const dimensionedScalar& rho1() const { return rho1_; }
        const dimensionedScalar& rho2() const { return rho2_; }
        //- Access Hertz-Knudsen gas constant [J/(kg K)]
        scalar gasConstant() const;
        //- Access evaporation accommodation coefficient []
        scalar evaporationCoeff() const;
        //- Access cached lattice heat capacity from the two-temperature model
        const dimensionedScalar& ClTTM() const { return ClTTM_; }

        //- Update cached lattice heat capacity from the two-temperature model
        void setClTTM(const dimensionedScalar& Cl);        
        //- Return const access to laser heating source
        const volScalarField& Q_laser() const
        {
            return Q_laser_;
        }
                //- Return access to laser heating source
        volScalarField& Q_laser()
        {
            return Q_laser_;
        }
        //- Copy external laser heating source into internal field
        void setQLaser(const volScalarField& src);
        //- Return const access to phase change source
        const volScalarField& phaseChangeSource() const
        {
            return phaseChangeSource_;
        }
        //- Return const access to the implicit relaxation coefficient
        const volScalarField& phaseChangeRelaxCoeff() const
        {
            return phaseChangeRelaxCoeff_;
        }
        //- Access dg/dt field [1/s]
        const volScalarField& dgdt() const
        {
            return dgdt_;
        }
        //- Access physical phase-change mass flux [kg/(m^2 s)]
        const volScalarField& phaseChangeMassFlux() const
        {
            return phaseChangeMassFlux_;
        }
        //- Return surface tension coefficient
        tmp<volScalarField> sigma() const;

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();
        
        //- Update mixture properties
        virtual void correct();
        
        //- Return the name of the thermo physics
        virtual word thermoName() const;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const;
        
        //- Compute phase change source term
        virtual void computePhaseChange();

        //- Compute and return mass-transfer rate dg/dt
        virtual tmp<volScalarField> computeMassTransfer() const;
 
        // Access to thermodynamic state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            //- Direct access to mixture enthalpy/internal energy is
            //  currently unsupported and returns a null field. A runtime
            //  warning is issued to discourage use.
            virtual volScalarField& he()
            {
                WarningInFunction
                    << "twoPhaseMixtureThermo::he() is unsupported." << endl;
                return const_cast<volScalarField&>(volScalarField::null());
            }

            //- Enthalpy/Internal energy [J/kg]
            //  Unsupported: returns a null field with a warning.
            virtual const volScalarField& he() const
            {
                WarningInFunction
                    << "twoPhaseMixtureThermo::he() const is unsupported." << endl;
                return volScalarField::null();
            }

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy for cell-set
            //  Uses a Newton iteration (up to 50 steps, tolerance 1e-6 on
            //  the maximum |ΔT| correction) seeded with T0 to invert the
            //  mixture enthalpy, emitting a warning if convergence is not
            //  achieved.
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            //  Uses a Newton iteration (up to 50 steps, tolerance 1e-6 on
            //  the maximum |ΔT| correction) seeded with T0 to invert the
            //  mixture enthalpy, emitting a warning if convergence is not
            //  achieved.
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;


        // Fields derived from thermodynamic state variables

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity using pressure and temperature
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Density from pressure and temperature
            virtual tmp<scalarField> rhoEoS
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const;

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;


        // Fields derived from transport state variables

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const;

            //- Thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const;


            //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;


    // IO

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
