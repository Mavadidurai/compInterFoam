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

// METAL PHASE PROPERTIES for LIFT Test Case
// Density is temperature-dependent using Boussinesq approximation for thermal expansion
// Representative solid/liquid titanium data (ρ ≈ 4500 kg/m3 at 300 K decreasing
// to ≈ 3400 kg/m3 in the high-temperature liquid state).
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
        Cp              650;            // J/kg·K - Heat capacity (liquid Ti, Keene 1993)
        Hf              0.0;            // J/kg - Heat of formation
        Tref            300;            // K - Reference temperature
        Href            0;              // J/kg - Reference enthalpy
    }
    
    transport
    {
        mu              2.35e-3;        // Pa·s - Dynamic viscosity (liquid Ti, Keene 1993)
        Pr              0.032;          // Prandtl number (consistent with Ti data)
    }

    equationOfState
    {
        rho0            4515;           // kg/m³ - Reference density at T0
        beta            7.6e-5;         // 1/K - Thermal expansion coefficient for Ti
        T0              300;            // K - Reference temperature
    }
}

// LIFT-specific properties
Tsol                1941.0;             // K - Solidus temperature
Tliq                1941.0;             // K - Liquidus temperature
Tvap                3560.0;             // K - Vaporisation temperature
hf    		    9.1e6;     	        // J/kg - Correct Ti vaporization latent heat (Palik)
                                        // Previous value 3.65e5 was fusion heat (25× too low)
kappa               17.2;               // W/m·K - Thermal conductivity

// ************************************************************************* //
