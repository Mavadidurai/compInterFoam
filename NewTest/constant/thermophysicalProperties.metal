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
// Density is now temperature dependent via a polynomial relation derived from
// representative solid/liquid titanium data (ρ ≈ 4500 kg/m3 at 300 K decreasing
// to ≈ 3400 kg/m3 in the high-temperature liquid state).
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rPolynomial;    // Polynomial density representation
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
        Cp              560;            // J/kg·K - Heat capacity
        Hf              0.0;            // J/kg - Heat of formation
        Tref            300;            // K - Reference temperature
        Href            0;              // J/kg - Reference enthalpy
    }
    
    transport
    {
        mu              2.25e-3;        // Pa·s - Dynamic viscosity
        Pr              0.03;           // Prandtl number (liquid metals)
    }
    
equationOfState
{
    // rPolynomial expects the specific volume fit 1/ρ = Σ Ci T^i (T in K)
    // Fitted to (T, ρ) = [(300 K, 4500 kg/m3), (1200 K, 4300 kg/m3),
    // (1941 K, 4100 kg/m3), (3200 K, 3400 kg/m3)]
    C (0.0002172169461 1.922024813e-08 -9.483330982e-12 3.433386724e-15 0);
}


}

// LIFT-specific properties
Tsol                1941.0;             // K - Solidus temperature
Tliq                1941.0;             // K - Liquidus temperature
Tvap                3560.0;             // K - Vaporisation temperature
hf                  3.65e5;              // J/kg - Latent heat of fusion
kappa               17.2;               // W/m·K - Thermal conductivity

// ************************************************************************* //
