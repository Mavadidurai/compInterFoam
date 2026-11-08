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

// Metal phase: titanium, treated as a liquid at constant density 4506 kg/m3. In an incompressible
// model the quantity that sets the launch is the mass per area of the film, rho*h. With the
// deposited thickness h = 71.4 nm this needs the room-temperature density of solid Ti, 4506 kg/m3.
// The measured liquid density (4193 kg/m3 at the melting point, Ozawa et al. 2017) at the same h
// would remove 7 % of the metal; the liquid's expansion (thicker film) is not represented.
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
        rho0            4506;           // kg/m3, solid Ti at 300 K (keeps the deposited mass per area)
        beta            0;              // constant metal density (verification choice): Boussinesq beta*dT ~ 0.4 at 5 kK is outside its validity and drives unphysical km/s expansion
        T0              300;            // Reference temperature [K]
    }
}


// ************************************************************************* //
