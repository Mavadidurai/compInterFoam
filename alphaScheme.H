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

ITstream compressionScheme
(
    compressionSchemes.found(alphaScheme)
  ? mesh.divScheme(divAlphaName)
  : ITstream
    (
        divAlphaName,
        tokenList
        {
            word("Gauss"),
            word("interfaceCompression"),
            alphaScheme,
            alphaControls.lookup<scalar>("cAlpha")
        }
    )
);
