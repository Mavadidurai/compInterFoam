// Helper for assembling implicit/explicit source terms in the alpha equation
// using the same constructs as the upstream compressibleInterFoam solver.

tmp<volScalarField::Internal> Sp
(
    volScalarField::Internal::New
    (
        "Sp",
        mesh,
        dimensionedScalar("zeroSp", dimless/dimTime, 0)
    )
);

tmp<volScalarField::Internal> Su
(
    volScalarField::Internal::New
    (
        "Su",
        mesh,
        dimensionedScalar("zeroSu", dimless/dimTime, 0)
    )
);

scalar sourceMaxMag = alphaControls.lookupOrDefault<scalar>("sourceMaxMag", GREAT);

scalarField& SpField = Sp.ref();
scalarField& SuField = Su.ref();

const scalarField& alpha1Field = alpha1.internalField();
const scalarField& dgdtField = dgdt.internalField();

forAll(SpField, celli)
{
    scalar rate = dgdtField[celli];

    if (sourceMaxMag < GREAT)
    {
        rate = Foam::max(Foam::min(rate, sourceMaxMag), -sourceMaxMag);
    }

    // dgdt > 0: evaporation (metal -> gas) removes alpha1 at rate dgdt;
    // dgdt < 0: condensation restores it.  Written in the bounded Weller
    // form used by upstream compressibleInterFoam: the implicit
    // coefficient must be NEGATIVE for MULES boundedness (the previous
    // positive "stabilizing" diagonal weakened the bound and injected a
    // spurious 2*Sp*alpha1 source).
    if (rate > 0.0)
    {
        // Sink proportional to alpha1: equals -rate at the current value
        // and shuts off smoothly as alpha1 -> 0.
        SpField[celli] -= rate/Foam::max(alpha1Field[celli], scalar(1e-4));
    }
    else if (rate < 0.0)
    {
        const scalar coeff =
            rate/Foam::max(1.0 - alpha1Field[celli], scalar(1e-4));

        SpField[celli] += coeff;   // negative (rate < 0)
        SuField[celli] -= coeff;   // positive source, bounded as alpha1 -> 1
    }
}

tmp<volScalarField> divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);
