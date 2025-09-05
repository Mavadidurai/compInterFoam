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
        type        symmetryPlane;
    }
    zMin
    {
        type        symmetryPlane;
    }
    xMin
    {
        type        symmetryPlane;
    }
    xMax
    {
        type        symmetryPlane;
    }
    yMin
    {
        type    constantAlphaContactAngle;
        theta0  2;           // degrees
        limit   gradient;
        value   uniform 1;   // metal wetting at the substrate side
    }
    yMax
    {
        type        inletOutlet;
        inletValue  uniform 0;
        value       uniform 0;
    }
}

// ************************************************************************* //
