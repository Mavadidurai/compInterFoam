extern Foam::Switch verbose;
verbose = runTime.controlDict().lookupOrDefault<Foam::Switch>("verbose", false);
bool LTS = fv::localEulerDdt::enabled(mesh);
tmp<volScalarField> trDeltaT;
if (LTS)
{
    if (verbose)
    {
        Info<< "Using LTS" << endl;
    }
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

if (verbose)
{
    Info<< "Reading field p_rgh\n" << endl;
}
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

if (verbose)
{
    Info<< "Reading field U\n" << endl;
}
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

if (verbose)
{
    Info<< "Reading/calculating face flux field phi\n" << endl;
}
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
if (verbose)
{
    Info<< "Reading field T\n" << endl;
}
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
if (verbose)
{
    Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
}
// Read transport properties
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        true   // ensure registration
    )
);
twoPhaseMixtureThermo mixture(U, phi);
// Set verbosity flag from controlDict now that runTime is available

{
    Foam::tmp<Foam::volScalarField> tsig = mixture.sigma();
    const Foam::volScalarField& sig = tsig();

    const Foam::dimensionedScalar sigmaMin = Foam::gMin(sig);
    const Foam::dimensionedScalar sigmaMax = Foam::gMax(sig);
    const Foam::dimensionedScalar sigmaAvg =
        Foam::fvc::domainIntegrate(sig)
       / Foam::dimensionedScalar
        (
            "totalVolume",
            Foam::dimVolume,
            Foam::gSum(mesh.V())
        );

    Foam::dimensionedScalar sigmaInput
    (
        "sigmaInput",
        Foam::dimForce/Foam::dimLength,
        0.0
    );

    if (mixture.found("sigma"))
    {
        if (mixture.isDict("sigma"))
        {
            mixture.subDict("sigma").lookup("sigma") >> sigmaInput;
        }
        else
        {
            mixture.lookup("sigma") >> sigmaInput;
        }
    }

    Info<< "surfaceTension input    = " << sigmaInput << nl
        << "surfaceTension [min,max,avg] = "
        << sigmaMin << ", " << sigmaMax << ", " << sigmaAvg << endl;
}

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

if (verbose)
{
    Info<< "Phase fractions - alpha.metal: min = " << min(alpha1).value()
        << ", max = " << max(alpha1).value() << endl;
    Info<< "Phase fractions - alpha.air: min = " << min(alpha2).value()
        << ", max = " << max(alpha2).value() << endl;

    // Verify phase naming consistency
    Info<< "Phase names - metal: " << mixture.phase1Name()
        << ", air: " << mixture.phase2Name() << endl;
}

if (verbose)
{
    Info<< "Reading thermophysical properties\n" << endl;
}

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

if (verbose)
{
    Info<< "Reading g\n" << endl;
}
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
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p_rgh, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
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

// Mass-transfer rate provided by twoPhaseMixtureThermo
const volScalarField& dgdt = mixture.dgdt();
if (verbose)
{
    Info<< "Mass-transfer rate - dgdt: min = " << min(dgdt).value()
        << ", max = " << max(dgdt).value() << endl;
}
IOobject alphaPhi10Header
(
    "alphaPhi10",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

const bool alphaRestart = alphaPhi10Header.typeHeaderOk<surfaceScalarField>();

if (alphaRestart && verbose)
{
    Info << "Restarting alpha flux field" << endl;
}

surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);
// Interface normal calculations
volVectorField n(fvc::grad(alpha1)/(mag(fvc::grad(alpha1)) + dimensionedScalar("small", dimless/dimLength, SMALL)));
surfaceVectorField nf(fvc::interpolate(n));
surfaceScalarField nHatf(nf & mesh.Sf());

// Create transportModel model
compressibleInterPhaseTransportModel transportModel
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

bool useAdvancedCapturing =
    runTime.controlDict().lookupOrDefault<bool>("useAdvancedInterfaceCapturing", true);
    
// Create LIFT-specific models
if (verbose)
{
    Info<< "Creating femtosecond laser model\n" << endl;
}

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
if (verbose)
{
    Info<< "Creating two-temperature model\n" << endl;
}
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
    twoTemperatureDict,
    alpha1
);
// Synchronise cached lattice heat capacity with the two-temperature model
mixture.setClTTM(ttm.Cl());

// Retrieve alpha solver controls from fvSolution
const dictionary& alphaSolverControls = mesh.solverDict("alpha.metal");

label nAlphaSubCycles
(
    alphaSolverControls.lookupOrDefault<label>("nAlphaSubCycles", 1)
);
label nAlphaCorr
(
    alphaSolverControls.lookupOrDefault<label>("nAlphaCorr", 2)
);
bool MULESCorr
(
    alphaSolverControls.lookupOrDefault<Foam::Switch>("MULESCorr", false)
);
bool alphaApplyPrevCorr
(
    alphaSolverControls.lookupOrDefault<Foam::Switch>("alphaApplyPrevCorr", false)
);
scalar icAlpha
(
    alphaSolverControls.lookupOrDefault<scalar>("icAlpha", 0.0)
);
scalar scAlpha
(
    alphaSolverControls.lookupOrDefault<scalar>("scAlpha", 0.0)
);
// Silence potential unused-variable warnings when the alpha equation
// is not executed (eg, post-processing builds)
(void)nAlphaSubCycles;
(void)nAlphaCorr;
(void)MULESCorr;
(void)alphaApplyPrevCorr;
(void)icAlpha;
(void)scAlpha;
// Suppress unused variable warnings for fvOptions
(void)fvOptions;

// Always create the interface-capturing model so a recoil-pressure field exists
autoPtr<advancedInterfaceCapturing> pInterfaceCapturing;
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

if (verbose)
{
    if (useAdvancedCapturing)
    {
        Info<< "Creating advanced interface capturing\n" << endl;
    }
    else
    {
        Info<< "Advanced interface capturing disabled; using basic recoil model\n" << endl;
    }
}

if (verbose)
{
    // ✅ CLEAR: Field validation with intuitive naming
    const dimensionedScalar totalMetalVOF = fvc::domainIntegrate(alpha1);
    const dimensionedScalar totalDomainVolume
    (
        "totalDomainVolume",
        dimVolume,
        gSum(mesh.V())
    );
    const scalar metalVolume_um3 = (totalMetalVOF*1e18).value();  // Convert m³ to μm³
    const scalar metalFraction = (totalMetalVOF/totalDomainVolume).value();

    Info<< "LIFT Field Validation:" << nl
        << "  Total domain volume: " << (totalDomainVolume*1e18).value() << " μm³" << nl
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
}
