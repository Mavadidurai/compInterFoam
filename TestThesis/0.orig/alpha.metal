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
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0";
    object      alpha.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

// SLM: Metal powder initially absent (set by setFieldsDict)
internalField   uniform 0;

boundaryField
{
    zMax
    {
        type        zeroGradient;
    }
    zMin
    {
        type        zeroGradient;
    }
    xMin
    {
        type            zeroGradient;
    }
    xMax
    {
        type            zeroGradient;
    }
    yMin
    {
        type            constantAlphaContactAngle;
        theta0          2;              // degrees (2–5 is fine)
        limit           gradient;       // or 'theta' if you prefer
        value           uniform 1;      // metal: 1 / air: 0
    }

    // TOP (atmosphere/open) — let α leave the domain
    yMax
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;      // metal: 0 / air: 1
    }
}

// ************************************************************************* //
