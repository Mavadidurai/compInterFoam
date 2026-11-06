extern Foam::Switch verbose;
verbose = runTime.controlDict().lookupOrDefault<Foam::Switch>("verbose", false);
bool LTS = fv::localEulerDdt::enabled(mesh);
tmp<volScalarField> trDeltaT;
tmp<volScalarField> trSubDeltaT;
if (LTS)
{
    // localEuler (LTS) discretizes every ddt() with a per-cell pseudo-time
    // step (fv::localEulerDdt looks up its own rDeltaT field and ignores
    // runTime.deltaT()), decoupling the fields from a single, physically
    // meaningful elapsed time. This solver gates laser pulses
    // (femtosecondLaserModel), phase-change kinetics, and two-temperature
    // sub-cycling on runTime.value()/deltaT, so running with LTS silently
    // desynchronizes what the fields actually evolved against from what
    // gets reported/used as physical time. LTS is only valid for
    // steady-state acceleration, never for this time-accurate solver.
    FatalErrorInFunction
        << "ddtSchemes selects localEuler (LTS), which is incompatible with "
        << "compInterFoam: laser pulse timing, phase-change kinetics, and "
        << "two-temperature sub-cycling all rely on a single physically "
        << "consistent runTime.deltaT() across cells. Use a time-accurate "
        << "scheme such as Euler or backward instead."
        << exit(FatalError);

    if (verbose && Foam::Pstream::master())
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
            dimensionedScalar("oneOverTime", dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    trSubDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rSubDeltaTName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("oneOverSubTime", dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    
}

if (verbose && Foam::Pstream::master())
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

if (verbose && Foam::Pstream::master())
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

if (verbose && Foam::Pstream::master())
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
/*
surfaceVectorField recoilTractionf
(
    IOobject
    (
        "recoilTractionf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "zeroRecoilTraction",
        dimForce/dimVolume,

        Foam::vector::zero
    )
);
*/
if (verbose && Foam::Pstream::master())
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

// Fix: compressibility dilatation source for the alpha equation [1/s],
// set in pEqn.H after the pressure solve (stock compressibleInterFoam dgdt).
volScalarField dgdtComp
(
    IOobject
    (
        "dgdtComp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dgdtComp", dimless/dimTime, 0.0)
);

volScalarField gasMetalHeatFlux
(
    IOobject
    (
        "gasMetalHeatFlux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zeroGasMetalHeatFlux",
        dimPower/dimVolume,
        0.0
    )
);
// Create the thermo system - this will create alpha.metal and alpha.air fields
if (verbose && Foam::Pstream::master())
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
(void)transportProperties;

word metalGroup(word::null);
word gasGroup(word::null);

if (transportProperties.found("phases"))
{
    wordList phases(transportProperties.lookup("phases"));

    if (phases.size() > 0)
    {
        metalGroup = phases[0];
    }

    if (phases.size() > 1)
    {
        gasGroup = phases[1];
    }
}

if (metalGroup.empty())
{
    metalGroup = "metal";
}

if (gasGroup.empty())
{
    gasGroup = "gas";
}

const word metalTName = IOobject::groupName("T", metalGroup);
const word gasTName = IOobject::groupName("T", gasGroup);

autoPtr<volScalarField> metalTPtr;
autoPtr<volScalarField> gasTPtr;

if (!mesh.foundObject<volScalarField>(metalTName))
{
    IOobject metalTHeader
    (
        metalTName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!metalTHeader.typeHeaderOk<volScalarField>())
    {
        metalTPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    metalTName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false   // Fix: not registered (see below)
                ),
                T
            )
        );

        metalTPtr->boundaryFieldRef() = T.boundaryField();
        metalTPtr->correctBoundaryConditions();
        metalTPtr->write();
        // Fix: the copy only seeds the file. Previously it stayed in the
        // registry, so on a fresh start the phase thermo found and kept this
        // frozen copy, while a restart read its own field from file. Freeing it
        // makes the thermo read the file and own the field, exactly as on a
        // restart (the path checked during the solver verification).
        metalTPtr.clear();
    }
}

if (!mesh.foundObject<volScalarField>(gasTName))
{
    IOobject gasTHeader
    (
        gasTName,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!gasTHeader.typeHeaderOk<volScalarField>())
    {
        gasTPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    gasTName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false   // Fix: not registered (see below)
                ),
                T
            )
        );

        gasTPtr->boundaryFieldRef() = T.boundaryField();
        gasTPtr->correctBoundaryConditions();
        gasTPtr->write();
        // Fix: the copy only seeds the file. Previously it stayed in the
        // registry, so on a fresh start the phase thermo found and kept this
        // frozen copy, while a restart read its own field from file. Freeing it
        // makes the thermo read the file and own the field, exactly as on a
        // restart (the path checked during the solver verification).
        gasTPtr.clear();
    }
}

twoPhaseMixtureThermo mixture(U, phi);
// Set verbosity flag from controlDict now that runTime is available
if (verbose && Foam::Pstream::master())
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

// Guard: the models look the phase fractions up by name
// (femtosecondLaserModel: "alpha.metal", advancedInterfaceCapturing:
// "alpha.air").  A duplicate registration ("already checked in") leaves the
// FIRST object in the registry; if that is not the mixture's own field, the
// models would read a stale copy.  Report which one the registry holds.
for (const volScalarField* a : {&alpha1, &alpha2})
{
    const volScalarField* reg = mesh.findObject<volScalarField>(a->name());
    if (reg == a)
    {
        Info<< "Phase-fraction registry check: registry '" << a->name()
            << "' is the mixture's phase fraction (OK)" << endl;
    }
    else
    {
        // fatal: the laser and interface models look these
        // names up and would silently compute sources from a stale field
        FatalErrorInFunction
            << "registry entry '" << a->name() << "' is "
            << (reg ? "a different (stale) object" : "missing")
            << ", not the mixture's phase fraction. Models that look it up"
            << " by name would read the wrong field." << exit(FatalError);
    }
}

const scalar metalFractionTol = 1e-6;
const scalar alphaMax = gMax(alpha1);

if (Foam::Pstream::master() && alphaMax < metalFractionTol)
{
    WarningInFunction
        << "alpha.metal is effectively empty (max(alpha.metal) = "
        << alphaMax
        << "). If this is an intentional vacuum shot, continue; otherwise"
        << " rerun setFields so the donor film matches the LIFT setup."
        << endl;
}

// Ensure proper field registration by always creating an old-time storage
// so that a cold start still has valid previous-time data for alpha
alpha1.oldTime();
alpha2.oldTime();

    if (verbose && Foam::Pstream::master())
{
    Info<< "Phase fractions - alpha.metal: min = " << min(alpha1).value()
        << ", max = " << max(alpha1).value() << endl;
    Info<< "Phase fractions - alpha.air: min = " << min(alpha2).value()
        << ", max = " << max(alpha2).value() << endl;

}

    if (verbose && Foam::Pstream::master())
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

    if (verbose && Foam::Pstream::master())
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

const Foam::vector gValue = g.value();
const Foam::scalar gMag = mag(gValue);

// ✅ FIXED: Proper dimensional calculation for ghRef preserving gravity direction
dimensionedScalar ghRef
(
    "ghRef",
    g.dimensions()*dimLength,
    gMag > SMALL
  ? (gValue & (-(gValue/gMag)*hRef.value()))
  : 0.0
);

volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);
// Fix: use the thermo's registered pressure (as stock compressibleInterFoam
// does) instead of a second, never-registered volScalarField named "p".
// Previously the solver's p = p_rgh + rho*gh failed to check in to the
// registry ("Attempt to checkIn object with name p which was already checked
// in"), so the thermo pressure -- used by he(p,T), the perfect-gas density and
// every lookupObject("p") -- stayed frozen at its initial value.
// p_rgh must therefore be an ABSOLUTE pressure (e.g. 1e5 Pa), not a gauge.
volScalarField& p = mixture.p();
p = p_rgh + rho*gh;
p.correctBoundaryConditions();
// Only the gas phase must have p > 0 (perfect-gas density); tension (p < 0)
// in the liquid metal is physical and is allowed (e.g. on restart).
// with an incompressible (rhoConst, psi = 0) gas, a local p < 0 on
// restart is only the incompressible pressure's response to the film pulling
// away from the gas (seen at 1.2x fluence); it has no thermodynamic effect,
// so it is a warning there. Fatal remains: an uninitialised/gauge p_rgh
// (gas-volume mean p <= 0), or p <= 0 anywhere with a compressible gas.
scalar minGasP = GREAT;
scalar sumGasPV = 0;
scalar sumGasV = 0;
forAll(p, celli)
{
    // positivity: every cell with a non-trivial gas fraction (the gas EOS is
    // also evaluated in mixed interface cells), as in the post-solve check in
    // pEqn.H; the gas-volume MEAN still uses majority-gas cells
    if (mixture.alpha2()[celli] > 1e-6)
    {
        minGasP = Foam::min(minGasP, p[celli]);
    }
    if (mixture.alpha2()[celli] > 0.5)
    {
        sumGasPV += p[celli]*mesh.V()[celli];
        sumGasV += mesh.V()[celli];
    }
}
reduce(minGasP, minOp<scalar>());
reduce(sumGasPV, sumOp<scalar>());
reduce(sumGasV, sumOp<scalar>());
const scalar meanGasP = sumGasPV/Foam::max(sumGasV, VSMALL);
const bool gasCompressible = gMax(mixture.thermo2().psi().primitiveField()) > 0;
// With no gas cell at all (single-phase verification cases) the gas mean is
// undefined: check the mean over the whole domain instead.
const scalar meanCheckP =
    sumGasV > VSMALL ? meanGasP : gSum((p.primitiveField()*mesh.V())())/gSum(mesh.V());
// on a restart (startTime > 0) a negative mean gas pressure is the
// incompressible gas's tension from the running solution, not an
// uninitialised gauge field, so it is reported as a warning instead.
const bool freshStart = (runTime.startTime().value() <= 0);
if ((freshStart && meanCheckP <= 0) || (sumGasV > VSMALL && minGasP <= 0 && gasCompressible))
{
    FatalErrorInFunction
        << "Absolute pressure p = p_rgh + rho*gh is <= 0 in the gas (min " << minGasP
        << " Pa, gas-volume mean " << meanGasP << " Pa). Since the pressure-field fix"
        << " p_rgh is an absolute pressure: initialise 0/p_rgh to the ambient"
        << " pressure (e.g. uniform 1e5) and set any totalPressure p0 /"
        << " pRefValue accordingly." << exit(FatalError);
}
else if (meanCheckP <= 0 || (sumGasV > VSMALL && minGasP <= 0))
{
    WarningInFunction
        << "Local absolute gas pressure <= 0 (min " << minGasP
        << " Pa, gas-volume mean " << meanGasP << " Pa) with an"
        << " incompressible gas: accepted (tension of the incompressible"
        << " pressure field, no thermodynamic effect)." << endl;
}
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p_rgh, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
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
    if (verbose && Foam::Pstream::master())
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

if (alphaRestart && verbose && Foam::Pstream::master())
{
    Info << "Restarting alpha flux field" << endl;
}

surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);
// Determine requested transport modelling approach and validate inputs
word transportSimulationType("laminar");

const auto readTransportSimulationType =
    [&](const word& dictName) -> bool
    {
        IOobject io
        (
            dictName,
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>())
        {
            return false;
        }

        IOdictionary dict(io);

        if (dict.found("simulationType"))
        {
            dict.lookup("simulationType") >> transportSimulationType;
            return true;
        }

        return false;
    };

if (!readTransportSimulationType("momentumTransport"))
{
    readTransportSimulationType("turbulenceProperties");
}

if (transportSimulationType == "twoPhaseTransport")
{
    const word phase1Name = alpha1.group();
    const word phase2Name = alpha2.group();

    const auto dictionaryExists =
        [&](const word& objectName) -> bool
        {
            IOobject io
            (
                objectName,
                runTime.constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );

            return io.typeHeaderOk<IOdictionary>();
        };

    const bool phase1DictFound =
        dictionaryExists(IOobject::groupName("momentumTransport", phase1Name))
     || dictionaryExists("phase1RASProperties")
     || dictionaryExists("phase1LESProperties")
     || dictionaryExists("phase1LaminarProperties");

    const bool phase2DictFound =
        dictionaryExists(IOobject::groupName("momentumTransport", phase2Name))
     || dictionaryExists("phase2RASProperties")
     || dictionaryExists("phase2LESProperties")
     || dictionaryExists("phase2LaminarProperties");

    if (!phase1DictFound || !phase2DictFound)
    {
        FatalErrorInFunction
            << "simulationType is set to twoPhaseTransport but the required phase"
            << " turbulence dictionaries are missing in '" << runTime.constant()
            << "'. Ensure that dictionaries such as 'momentumTransport." << phase1Name
            << "' and 'momentumTransport." << phase2Name
            << "' (or the corresponding phase1*/phase2* RAS/LES properties) are"
            << " provided." << nl
            << exit(FatalError);
    }
}

// the interface normal is taken from mixture.nHatf() (updated at
// every alpha correction). A static copy built here from the initial alpha
// was used by alphaEqn.H and froze the compression at the initial interface.

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

Switch useAdvancedCapturing
(
    runTime.controlDict().lookupOrDefault<Switch>
    (
        "useAdvancedInterfaceCapturing",
        true
    )
);
    
// Create LIFT-specific models
if (verbose && Foam::Pstream::master())
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
const dictionary& globalControlDict = runTime.controlDict();

const bool hasLaserStartControl = globalControlDict.found("laserStartTime");
const bool hasLaserEndControl = globalControlDict.found("laserEndTime");
const bool hasLaserEndDict = laserPropertiesDict.found("laserEndTime");

const scalar laserStartTimeControl =
    hasLaserStartControl
  ? readScalar(globalControlDict.lookup("laserStartTime"))
  : laserPropertiesDict.lookupOrDefault<scalar>("laserStartTime", 0.0);
laserPropertiesDict.set("laserStartTime", laserStartTimeControl);

scalar laserEndTimeDefault =
    hasLaserEndDict
  ? readScalar(laserPropertiesDict.lookup("laserEndTime"))
  : laserStartTimeControl;

const scalar laserEndTimeControl =
    hasLaserEndControl
  ? readScalar(globalControlDict.lookup("laserEndTime"))
  : laserEndTimeDefault;
laserPropertiesDict.set("laserEndTime", laserEndTimeControl);

if (!hasLaserEndControl && !hasLaserEndDict)
{
    FatalIOErrorInFunction(globalControlDict)
        << "Specify 'laserEndTime' in controlDict or constant/laserProperties" << nl
        << exit(FatalIOError);
}
femtosecondLaserModel laser(mesh, laserPropertiesDict);
if (verbose && Foam::Pstream::master())
{
    Info<< "Creating two-temperature model\n" << endl;
}
// README documents twoTemperatureProperties as readable from either
// system/controlDict or the optional constant/twoTemperatureProperties;
// the latter was never actually wired up. Read constant/ as the base (if
// present) and let a controlDict.twoTemperatureProperties subdict override
// it, matching how laserProperties/laserEndTime are merged above.
IOobject twoTemperaturePropertiesHeader
(
    "twoTemperatureProperties",
    runTime.constant(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
);

const bool haveTwoTemperaturePropertiesFile =
    twoTemperaturePropertiesHeader.typeHeaderOk<IOdictionary>();

dictionary twoTemperatureDict;

if (haveTwoTemperaturePropertiesFile)
{
    twoTemperatureDict = IOdictionary(twoTemperaturePropertiesHeader);
}

if (runTime.controlDict().found("twoTemperatureProperties"))
{
    twoTemperatureDict.merge
    (
        runTime.controlDict().subDict("twoTemperatureProperties")
    );
}
else if (!haveTwoTemperaturePropertiesFile)
{
    FatalIOErrorInFunction(runTime.controlDict())
        << "Specify 'twoTemperatureProperties' in controlDict or"
        << " constant/twoTemperatureProperties"
        << exit(FatalIOError);
}

scalar energyToleranceValue =
    twoTemperatureDict.lookupOrDefault<scalar>("energyTolerance", 0.1);
energyToleranceValue =
    twoTemperatureDict.lookupOrDefault<scalar>("energyTol", energyToleranceValue);
twoTemperatureDict.set("energyTolerance", energyToleranceValue);

dimensionedScalar energyTolerance
(
    "energyTolerance",
    dimless,
    energyToleranceValue
);

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
label nLimiterIter
(
    alphaSolverControls.lookupOrDefault<label>("nLimiterIter", 3)
);
bool MULESCorr
(
    alphaSolverControls.lookupOrDefault<Foam::Switch>("MULESCorr", false)
);
bool alphaApplyPrevCorr
(
    alphaSolverControls.lookupOrDefault<Foam::Switch>("alphaApplyPrevCorr", false)
);
// Persists across time steps (alphaEqn.H is #included afresh each call):
// caches the previous MULESCorr compression flux so it can seed the next
// correction when alphaApplyPrevCorr is enabled. Left invalid until the
// first MULESCorr pass populates it.
tmp<surfaceScalarField> talphaPhi1Corr0;
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
(void)nLimiterIter;
(void)MULESCorr;
(void)alphaApplyPrevCorr;
(void)icAlpha;
(void)scAlpha;
// Suppress unused variable warnings for fvOptions
(void)fvOptions;

// Conditionally create the interface-capturing model so legacy behaviour is preserved
autoPtr<advancedInterfaceCapturing> pInterfaceCapturing;
autoPtr<volScalarField> legacyRecoilPressure;
const volScalarField* recoilPressurePtr = nullptr;
if (useAdvancedCapturing)
{
    if (verbose && Foam::Pstream::master())

    {
        Info<< "Creating advanced interface capturing\n" << endl;
    }

    pInterfaceCapturing.reset
    (
        new advancedInterfaceCapturing
        (
            mesh,
            alpha1,
            mixture,
            T,
            &ttm.Tl()
        )
    );
    recoilPressurePtr = &pInterfaceCapturing->recoilPressure();
}
else
{
    if (verbose && Foam::Pstream::master())
    {
        Info<< "Advanced interface capturing disabled; using basic recoil model\n"
            << endl;
    }

    legacyRecoilPressure.reset
    (
        new volScalarField
        (
            IOobject
            (
                "recoilPressure",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    );

    legacyRecoilPressure->correctBoundaryConditions();
    // get() keeps ownership with the autoPtr (ptr() would release it and
    // leave legacyRecoilPressure invalid for the refresh/write calls).
    recoilPressurePtr = legacyRecoilPressure.get();
}
autoPtr<enhancedLIFTPhysics> liftPhysics;

if (runTime.controlDict().found("enhancedLIFTPhysics"))
{
    if (verbose && Foam::Pstream::master())
    {
        Info<< "Creating enhanced LIFT physics models\n" << endl;
    }

    liftPhysics.reset
    (
        new enhancedLIFTPhysics
        (
            mesh,
            mixture,
            runTime.controlDict().subDict("enhancedLIFTPhysics")
        )
    );
}
(void)recoilPressurePtr;

if (verbose && Foam::Pstream::master())
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
