// CORRECTED CODE FOR compInterFoam.C (OpenFOAM v2406)
// Replace lines 223-264 with this version
// Uses modern get<> API - NO deprecation warnings

// ===== Extract phase change parameters =====
if (controlDict.found("phaseChangeCoeffs"))
{
    const Foam::dictionary& phaseChangeDict = controlDict.subDict("phaseChangeCoeffs");

    // Dimensioned parameters - use modern get<> API (OpenFOAM v2406)
    if (phaseChangeDict.found("Tvap"))
    {
        input_Tvap = phaseChangeDict.get<Foam::dimensionedScalar>("Tvap").value();
    }
    if (phaseChangeDict.found("Tsol"))
    {
        input_Tsol = phaseChangeDict.get<Foam::dimensionedScalar>("Tsol").value();
    }
    if (phaseChangeDict.found("hf"))
    {
        input_hf = phaseChangeDict.get<Foam::dimensionedScalar>("hf").value();
    }
    if (phaseChangeDict.found("gasConstant"))
    {
        input_gasConstant = phaseChangeDict.get<Foam::dimensionedScalar>("gasConstant").value();
    }
    if (phaseChangeDict.found("p_ref"))
    {
        input_p_ref = phaseChangeDict.get<Foam::dimensionedScalar>("p_ref").value();
    }
    if (phaseChangeDict.found("maxSource"))
    {
        input_maxSource = phaseChangeDict.get<Foam::dimensionedScalar>("maxSource").value();
    }

    // Dimensionless parameters - read directly as scalar
    input_evaporationCoeff = phaseChangeDict.lookupOrDefault<Foam::scalar>("evaporationCoeff", 0.0);
    input_relaxationTime = phaseChangeDict.lookupOrDefault<Foam::scalar>("relaxationTime", 1e-11);
    input_alphaMin = phaseChangeDict.lookupOrDefault<Foam::scalar>("alphaMin", 0.001);
    input_metalFractionCutoff = phaseChangeDict.lookupOrDefault<Foam::scalar>("metalFractionCutoff", 1e-6);
}

// ===== Extract two-temperature model parameters =====
if (controlDict.found("twoTemperatureProperties"))
{
    const Foam::dictionary& twoTempDict = controlDict.subDict("twoTemperatureProperties");

    // Dimensioned parameters - use modern get<> API (OpenFOAM v2406)
    if (twoTempDict.found("Cl"))
    {
        input_Cl = twoTempDict.get<Foam::dimensionedScalar>("Cl").value();
    }
    if (twoTempDict.found("De"))
    {
        input_De = twoTempDict.get<Foam::dimensionedScalar>("De").value();
    }

    // Extract Ce coefficient from polynomial if present
    if (twoTempDict.found("Ce"))
    {
        const Foam::dictionary& CeDict = twoTempDict.subDict("Ce");
        if (CeDict.found("coeffs"))
        {
            Foam::List<Foam::Tuple2<Foam::scalar, Foam::scalar>> coeffs =
                CeDict.lookup("coeffs");
            // Look for linear coefficient (second entry)
            if (coeffs.size() > 1)
            {
                input_Ce_coeff = coeffs[1].first();
            }
        }
    }

    // Dimensionless parameters - read directly as scalar
    input_maxTe = twoTempDict.lookupOrDefault<Foam::scalar>("maxTe", 20000.0);
    input_maxTl = twoTempDict.lookupOrDefault<Foam::scalar>("maxTl", 10000.0);
    input_minTe = twoTempDict.lookupOrDefault<Foam::scalar>("minTe", 200.0);
}

// Optional: Add validation output after extraction
if (Foam::Pstream::master())
{
    Info<< "CSV Parameter Extraction Validation:" << Foam::nl
        << "  Tvap = " << input_Tvap << " K" << Foam::nl
        << "  Tsol = " << input_Tsol << " K" << Foam::nl
        << "  hf = " << input_hf << " J/kg" << Foam::nl
        << "  gasConstant = " << input_gasConstant << " J/(kg·K)" << Foam::nl
        << "  p_ref = " << input_p_ref << " Pa" << Foam::nl
        << "  maxSource = " << input_maxSource << " W/m³" << Foam::nl
        << "  Cl = " << input_Cl << " J/(m³·K)" << Foam::nl
        << "  De = " << input_De << " m²/s" << Foam::endl;
}
