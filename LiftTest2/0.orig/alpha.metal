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
      back 		{ 	 type        symmetryPlane;	}
      front   		{        type        symmetryPlane;     }
      left  		{        type        symmetryPlane;     }
      right    		{        type        symmetryPlane;     }
      donorSubstrate     {	 type        zeroGradient;      }  // ✅ CORRECTED
    receiver
    {
        type        zeroGradient;  // ✅ CORRECTED
    }
}

// ************************************************************************* //
