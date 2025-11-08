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
const scalarField& alpha2Field = alpha2.internalField();
const scalarField& rho1Field = rho1.internalField();
const scalarField& rho2Field = rho2.internalField();
const scalarField& dgdtField = dgdt.internalField();

forAll(SpField, celli)
{
    scalar rate = dgdtField[celli];

    if (sourceMaxMag < GREAT)
    {
        rate = Foam::max(Foam::min(rate, sourceMaxMag), -sourceMaxMag);
    }

    const scalar a1 = alpha1Field[celli];
    const scalar a2 = alpha2Field[celli];

    const scalar rhoMix = Foam::max
    (
        a1*rho1Field[celli] + a2*rho2Field[celli],
        SMALL
    );

    const scalar absRate = Foam::mag(rate);
    const scalar evapRate = 0.5*(rate + absRate);
    const scalar condRate = 0.5*(absRate - rate);

    const scalar evapWeight = (a2*rho2Field[celli])/rhoMix;
    const scalar condWeight = (a1*rho1Field[celli])/rhoMix;

    const scalar diag = -(evapRate*evapWeight + condRate*condWeight);

    SpField[celli] = diag;
    SuField[celli] = -rate - diag*a1;
}

tmp<volScalarField> divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);
