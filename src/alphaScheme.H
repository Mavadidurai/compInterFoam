// This header mirrors the upstream OpenFOAM alphaScheme helper so that the
// local alpha equation follows the standard structure and can pick up future
// updates with minimal churn.

static const wordHashSet compressionSchemes
{
    "interfaceCompression",
    "noInterfaceCompression",
    "PLIC",
    "PLICU",
    "MPLIC",
    "MPLICU"
};

static const word divAlphaName("div(phi,alpha)");

const word alphaScheme(mesh.divScheme(divAlphaName)[1].wordToken());

tokenList defaultCompressionScheme(4);
defaultCompressionScheme[0] = token(word("Gauss"));
defaultCompressionScheme[1] = token(word("interfaceCompression"));
defaultCompressionScheme[2] = token(alphaScheme);
defaultCompressionScheme[3] = Foam::token
(
    alphaControls.lookupOrDefault<scalar>("cAlpha", mixture.cAlpha())
);

ITstream compressionScheme
(
    compressionSchemes.found(alphaScheme)
  ? mesh.divScheme(divAlphaName)
  : ITstream
    (
        divAlphaName,
        defaultCompressionScheme
    )
);
word compressionFluxSchemeName(alphaScheme);

{
    word integrationScheme(alphaScheme);
    compressionScheme >> integrationScheme;

    if (compressionScheme.good())
    {
        compressionScheme >> compressionFluxSchemeName;
    }
    else
    {
        compressionFluxSchemeName = integrationScheme;
    }
}

if (!compressionSchemes.found(compressionFluxSchemeName))
{
    compressionFluxSchemeName = word(defaultCompressionScheme[1].wordToken());
}

compressionScheme.rewind();
