"""Read material values from a case, so the analysis tools do not hard-code them."""
import os, re


def metal_density(case, default=None):
    """rho0 of the metal phase from constant/thermophysicalProperties.metal [kg/m3]."""
    p = os.path.join(case, "constant", "thermophysicalProperties.metal")
    try:
        s = re.sub(r"//.*", "", open(p).read())
        m = re.search(r"\brho0\s+([0-9.eE+-]+)\s*;", s)
        if m:
            return float(m.group(1))
    except OSError:
        pass
    if default is None:
        raise SystemExit("cannot read rho0 from %s" % p)
    return default


def surface_tension(case, default=None):
    """sigma from constant/thermophysicalProperties [N/m]."""
    p = os.path.join(case, "constant", "thermophysicalProperties")
    try:
        s = re.sub(r"//.*", "", open(p).read())
        m = re.search(r"^\s*sigma\s+(?:sigma\s+\[[^\]]*\]\s+)?([0-9.eE+-]+)\s*;", s, re.M)
        if m:
            return float(m.group(1))
    except OSError:
        pass
    if default is None:
        raise SystemExit("cannot read sigma from %s" % p)
    return default
