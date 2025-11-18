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
// CRITICAL FIX: Changed from Boussinesq to perfectGas to allow density collapse
// during superheat → vapor → plume. Boussinesq was preventing expansion needed for jetting.
// NOTE: This is a numerical approximation for compressibility, not physically perfect for Ti.
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState perfectGas;
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
        Cp              650;            // J/kg·K - Heat capacity (molten Ti, Keene 1993)
        Hf              0.0;            // J/kg - Heat of formation
        Tref            300;            // K - Reference temperature
        Href            0;              // J/kg - Reference enthalpy
    }

    transport
    {
        mu              2.35e-3;        // Pa·s - Dynamic viscosity (Keene 1993)
        Pr              0.032;          // Prandtl number (consistent with Ti data)
    }
    
    equationOfState
    {
        // perfectGas EOS parameters (numerical approximation for Ti compressibility)
        // rho0 kept for reference; solver uses psi/thermo for dynamic density
        rho0            4515;           // kg/m3 at reference (informational only)
    }
}

// LIFT-specific properties
Tsol                1941.0;             // K - Solidus temperature
Tliq                1941.0;             // K - Liquidus temperature
Tvap                2200.0;             // K - Superheated liquid regime (CORRECTED from 3560)
hf                  9.1e6;             // Latent heat of vaporisation for Ti (J/kg)
kappa               0;                  // CORRECTED: Let TTM control thermal conduction via ke and kl

// ************************************************************************* //
