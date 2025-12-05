// CORRECTED CODE FOR compInterFoam.C
// Replace lines 223-264 with this corrected version

// ===== Extract phase change parameters =====
if (controlDict.found("phaseChangeCoeffs"))
{
    const Foam::dictionary& phaseChangeDict = controlDict.subDict("phaseChangeCoeffs");

    // Dimensioned parameters - must extract from dimensionedScalar objects
    if (phaseChangeDict.found("Tvap"))
    {
        try
        {
            Foam::dimensionedScalar Tvap_dim(phaseChangeDict.lookup("Tvap"));
            input_Tvap = Tvap_dim.value();
        }
        catch (...)
        {
            input_Tvap = phaseChangeDict.lookupOrDefault<Foam::scalar>("Tvap", 0.0);
        }
    }

    if (phaseChangeDict.found("Tsol"))
    {
        try
        {
            Foam::dimensionedScalar Tsol_dim(phaseChangeDict.lookup("Tsol"));
            input_Tsol = Tsol_dim.value();
        }
        catch (...)
        {
            input_Tsol = phaseChangeDict.lookupOrDefault<Foam::scalar>("Tsol", 0.0);
        }
    }

    if (phaseChangeDict.found("hf"))
    {
        try
        {
            Foam::dimensionedScalar hf_dim(phaseChangeDict.lookup("hf"));
            input_hf = hf_dim.value();
        }
        catch (...)
        {
            input_hf = phaseChangeDict.lookupOrDefault<Foam::scalar>("hf", 0.0);
        }
    }

    if (phaseChangeDict.found("gasConstant"))
    {
        try
        {
            Foam::dimensionedScalar gasConstant_dim(phaseChangeDict.lookup("gasConstant"));
            input_gasConstant = gasConstant_dim.value();
        }
        catch (...)
        {
            input_gasConstant = phaseChangeDict.lookupOrDefault<Foam::scalar>("gasConstant", 0.0);
        }
    }

    if (phaseChangeDict.found("p_ref"))
    {
        try
        {
            Foam::dimensionedScalar p_ref_dim(phaseChangeDict.lookup("p_ref"));
            input_p_ref = p_ref_dim.value();
        }
        catch (...)
        {
            input_p_ref = phaseChangeDict.lookupOrDefault<Foam::scalar>("p_ref", 101325.0);
        }
    }

    if (phaseChangeDict.found("maxSource"))
    {
        try
        {
            Foam::dimensionedScalar maxSource_dim(phaseChangeDict.lookup("maxSource"));
            input_maxSource = maxSource_dim.value();
        }
        catch (...)
        {
            input_maxSource = phaseChangeDict.lookupOrDefault<Foam::scalar>("maxSource", 0.0);
        }
    }

    // Dimensionless parameters - read directly as scalar (these are correct)
    input_evaporationCoeff = phaseChangeDict.lookupOrDefault<Foam::scalar>("evaporationCoeff", 0.0);
    input_relaxationTime = phaseChangeDict.lookupOrDefault<Foam::scalar>("relaxationTime", 1e-11);
    input_alphaMin = phaseChangeDict.lookupOrDefault<Foam::scalar>("alphaMin", 0.001);
    input_metalFractionCutoff = phaseChangeDict.lookupOrDefault<Foam::scalar>("metalFractionCutoff", 1e-6);
}

// ===== Extract two-temperature model parameters =====
if (controlDict.found("twoTemperatureProperties"))
{
    const Foam::dictionary& twoTempDict = controlDict.subDict("twoTemperatureProperties");

    // Dimensioned parameters - must extract from dimensionedScalar objects
    if (twoTempDict.found("Cl"))
    {
        try
        {
            Foam::dimensionedScalar Cl_dim(twoTempDict.lookup("Cl"));
            input_Cl = Cl_dim.value();
        }
        catch (...)
        {
            input_Cl = twoTempDict.lookupOrDefault<Foam::scalar>("Cl", 0.0);
        }
    }

    if (twoTempDict.found("De"))
    {
        try
        {
            Foam::dimensionedScalar De_dim(twoTempDict.lookup("De"));
            input_De = De_dim.value();
        }
        catch (...)
        {
            input_De = twoTempDict.lookupOrDefault<Foam::scalar>("De", 0.0);
        }
    }

    // Extract Ce coefficient from polynomial if present (this is already correct)
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

    // Dimensionless parameters (these are correct)
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
