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
defaultCompressionScheme[3] = token(alphaControls.lookup<scalar>("cAlpha"));

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
