/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Metal phase: titanium, treated as a liquid at constant density 4515 kg/m3
// (beta = 0 below switches the thermal expansion off).
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState Boussinesq;
    specie          specie;
    energy          sensibleEnthalpy;

}

mixture
{
    specie
    {
        molWeight       47.867;         // kg/kmol - Titanium molecular weight
    }

    thermodynamics
   {
        Cp              650;            // J/kg·K - Heat capacity
        Hf              0.0;            // J/kg - Heat of formation
        Tref            300;            // K - Reference temperature
        Href            0;              // J/kg - Reference enthalpy
    }

    transport
    {
        mu              2.35e-3;        // Pa·s - Dynamic viscosity
        Pr              0.032;           // Prandtl number (liquid metals)
    }

    equationOfState
    {
        rho0            4515;           // kg/m3 at reference temperature
        beta            0;              // constant metal density (verification choice): Boussinesq beta*dT ~ 0.4 at 5 kK is outside its validity and drives unphysical km/s expansion
        T0              300;            // Reference temperature [K]
    }
}

// LIFT-specific properties
Tsol                1941.0;             // K - Solidus temperature
Tliq                1941.0;             // K - Liquidus temperature
Tvap                3560.0;             // K - Vaporisation temperature
hf                  9.1e6;             // Latent heat of vaporisation for Ti (J/kg)
kappa               17.2;               // W/m·K - Thermal conductivity

// ************************************************************************* //
