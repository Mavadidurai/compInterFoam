
bool LTS = fv::localEulerDdt::enabled(mesh);
tmp<volScalarField> trDeltaT;
if (LTS)
{
    Info<< "Using LTS" << endl;
    trDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
}

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);
Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("T", dimTemperature, 300.0)
);
// Create the thermo system - this will create alpha.metal and alpha.air fields
Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

// ✅ INTUITIVE: Access alpha fields using clear naming
// The mixture.alpha1() refers to the first phase (metal)
// The mixture.alpha2() refers to the second phase (air)
volScalarField& alpha1 = mixture.alpha1();  // This is alpha.metal
volScalarField& alpha2 = mixture.alpha2();  // This is alpha.air

// Ensure proper field registration
if (alpha1.timeIndex() > 0)
{
    alpha1.oldTime();
    alpha2.oldTime();
}

Info<< "Phase fractions - alpha.metal: min = " << min(alpha1).value() 
    << ", max = " << max(alpha1).value() << endl;
Info<< "Phase fractions - alpha.air: min = " << min(alpha2).value() 
    << ", max = " << max(alpha2).value() << endl;

// Verify phase naming consistency
Info<< "Phase names - metal: " << mixture.phase1Name() 
    << ", air: " << mixture.phase2Name() << endl;

Info<< "Reading thermophysical properties\n" << endl;

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

Info<< "Reading g\n" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar hRef
(
    "hRef",
    dimLength,
    mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("hRef", scalar(0))
);

// ✅ FIXED: Proper dimensional calculation for ghRef
dimensionedScalar ghRef
(
    "ghRef",
    g.dimensions()*dimLength,
    mag(g.value()) > SMALL
  ? mag(g.value()) * hRef.value()
  : 0.0
);

volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);
surfaceScalarField rho1f(fvc::interpolate(rho1));
surfaceScalarField rho2f(fvc::interpolate(rho2));

surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt
(
    IOobject
    (
        "dgdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dgdt", dimless/dimTime, 0)
);

IOobject alphaPhi10Header
(
    "alphaPhi10",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

const bool alphaRestart = alphaPhi10Header.typeHeaderOk<surfaceScalarField>();

if (alphaRestart)
{
    Info << "Restarting alpha flux field" << endl;
}

surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);

tmp<surfaceScalarField> talphaPhi1Corr0;

// Alpha flux field
surfaceScalarField alphaPhi1
(
    IOobject
    (
        "alphaPhi1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

// Interface normal calculations
volVectorField n(fvc::grad(alpha1)/(mag(fvc::grad(alpha1)) + dimensionedScalar("small", dimless/dimLength, SMALL)));
surfaceVectorField nf(fvc::interpolate(n));
surfaceScalarField nHatf(nf & mesh.Sf());

// Create turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi10,
    mixture
);

volScalarField K("K", 0.5*magSqr(U));

IOMRFZoneList MRF(mesh);

fv::options& fvOptions(fv::options::New(mesh));

// Read transport properties
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

bool useAdvancedCapturing = 
    runTime.controlDict().lookupOrDefault<bool>("useAdvancedInterfaceCapturing", false);
    
// Create LIFT-specific models
Info<< "Creating femtosecond laser model\n" << endl;
IOdictionary laserPropertiesDict
(
    IOobject
    (
        "laserProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);
femtosecondLaserModel laser(mesh, laserPropertiesDict);

Info<< "Creating two-temperature model\n" << endl;
dictionary twoTemperatureDict
(
    runTime.controlDict().subDict("twoTemperatureProperties")
);

scalar energyTolerance =
    twoTemperatureDict.lookupOrDefault<scalar>("energyTolerance", 0.1);
energyTolerance =
    twoTemperatureDict.lookupOrDefault<scalar>("energyTol", energyTolerance);
twoTemperatureDict.set("energyTolerance", energyTolerance);

twoTemperatureModel ttm
(
    mesh,
    twoTemperatureDict
);


label nAlphaSubCycles(1);
label nAlphaCorr(2);
bool MULESCorr(false);
bool alphaApplyPrevCorr(false);
scalar icAlpha(0.0);
scalar scAlpha(0.0);

// Suppress unused variable warnings
(void)nAlphaSubCycles;
(void)nAlphaCorr;
(void)MULESCorr;
(void)alphaApplyPrevCorr;
(void)icAlpha;
(void)scAlpha;
(void)fvOptions;

// Pressure reference
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p_rgh, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
mesh.setFluxRequired(p_rgh.name());

// Create advanced interface capturing if needed
autoPtr<advancedInterfaceCapturing> pInterfaceCapturing;

if (useAdvancedCapturing)
{
    Info<< "Creating advanced interface capturing\n" << endl;
    pInterfaceCapturing.reset
    (
        new advancedInterfaceCapturing
        (
            mesh,
            alpha1,
            phi,
            mixture,
            T
        )
    );
}
else
{
    Info<< "Using standard interface capturing\n" << endl;
}

// ✅ CLEAR: Field validation with intuitive naming
dimensionedScalar totalMetalVOF = fvc::domainIntegrate(alpha1);
scalar totalDomainVolume = gSum(mesh.V());
scalar metalVolume_um3 = totalMetalVOF.value() * 1e18;  // Convert m³ to μm³
scalar metalFraction = totalMetalVOF.value() / totalDomainVolume;

Info<< "LIFT Field Validation:" << nl
    << "  Total domain volume: " << totalDomainVolume*1e18 << " μm³" << nl
    << "  Total metal volume: " << metalVolume_um3 << " μm³" << nl  
    << "  Metal fraction: " << metalFraction*100 << "%" << endl;

if (metalFraction > 0.001)  // If some metal detected
{
    Info<< "✅ Metal film detected: " << metalFraction*100 << "% of domain" << endl;
}
else
{
    WarningInFunction
        << "Very low metal content detected!" << nl
        << "  Found: " << metalFraction*100 << "% (" << metalVolume_um3 << " μm³)" << nl
        << "  Make sure to run 'setFields' to initialize the metal film!" << endl;
}

Info<< "Field setup completed successfully!" << endl;
