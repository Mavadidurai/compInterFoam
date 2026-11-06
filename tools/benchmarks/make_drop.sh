#!/bin/bash
# Build an axisymmetric drop benchmark case from the reference case (cases/fsLIFT_332mJ_axisymmetric)
# with the laser, phase change, recoil, vapour pocket and gravity switched off. Only the interface
# mechanics remain: surface tension, inertia and viscosity of a liquid titanium drop in gas.
#   make_drop.sh <dir> <cell size nm> <b/a aspect (1 = sphere)> <endTime s> <maxDeltaT s>
# Drop: Ti melt (rho 4515, mu 2.35e-3, sigma 1.64), R = 0.2 um, in rhoConst air, 300 K, g = 0.
set -e
D=$1; DX=$2; ASP=$3; TEND=$4; DTMAX=$5
SRC=$(cd "$(dirname "$0")/../../cases/fsLIFT_332mJ_axisymmetric" && pwd)
R=0.2; L=0.6
rm -rf "${D:?}"; mkdir -p $D/system
cp -r $SRC/0.orig $SRC/constant $D/; rm -rf $D/constant/polyMesh
for f in controlDict fvSchemes fvSolution decomposeParDict; do cp $SRC/system/$f $D/system/; done
N=$(python3 -c "print(round($L*1000/$DX))")
python3 - "$D" "$L" "$N" <<'PY'
import sys, math
D,L,N=sys.argv[1],float(sys.argv[2]),int(sys.argv[3])
th=math.radians(2.5); c,s=math.cos(th),math.sin(th)
x=L*c; z=L*s
V=f"""FoamFile {{ version 2.0; format ascii; class dictionary; object blockMeshDict; }}
scale 1e-6;
vertices
(
    (0 {-L} 0) ({x} {-L} {z}) ({x} {-L} {-z})
    (0 {L} 0)  ({x} {L} {z})  ({x} {L} {-z})
);
blocks ( hex (0 1 2 0 3 4 5 3) ({N} 1 {2*N}) simpleGrading (1 1 1) );
edges ();
boundary
(
    front     {{ type wedge; faces ((0 3 4 1)); }}
    back      {{ type wedge; faces ((0 2 5 3)); }}
    right     {{ type patch; faces ((1 4 5 2)); }}
    axis      {{ type empty; faces ((0 3 3 0)); }}
    substrate {{ type wall;  faces ((0 1 2 0)); }}
    donor     {{ type wall;  faces ((3 5 4 3)); }}
);
"""
open(D+'/system/blockMeshDict','w').write(V)
PY
# initial interface: smooth (tanh over ~1 cell) ellipsoid with semi-axes a (r), b = ASP*a (y), volume = sphere of R
python3 - "$D" "$R" "$ASP" "$DX" <<'PY'
import sys
D,R,asp,dx=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])
a=R*1e-6/asp**(1/3); b=asp*a; w=dx*1e-9
E=f"sqrt(pow(pos().x()/{a:.6e},2)+pow(pos().y()/{b:.6e},2)+pow(pos().z()/{a:.6e},2))"
S=f"""FoamFile {{ version 2.0; format ascii; class dictionary; object setExprFieldsDict; }}
expressions
(
    alphaMetal {{ field alpha.metal; dimensions [0 0 0 0 0 0 0]; expression #{{ 0.5*(1 - tanh(({E}-1)*{a:.6e}/{w:.6e})) #}}; }}
    alphaAir   {{ field alpha.air;   dimensions [0 0 0 0 0 0 0]; expression #{{ 0.5*(1 + tanh(({E}-1)*{a:.6e}/{w:.6e})) #}}; }}
);
"""
open(D+'/system/setExprFieldsDict','w').write(S)
PY
cd $D
printf 'FoamFile { version 2.0; format ascii; class dictionary; location "constant"; object fvOptions; }\n// benchmark: substrate constraints removed (no substrate zone in the drop mesh)\n' > constant/fvOptions
sed -i 's/^value .*/value           (0 0 0);   \/\/ benchmark: no gravity/' constant/g
sed -i 's/^pulseEnergy .*/pulseEnergy              1e-30;   \/\/ benchmark: laser off/' constant/laserProperties
sed -i -E 's/^(\s*tEnd\s*)\(5e-8\);.*/\1(1e-30);   \/\/ benchmark: phase change off/' system/controlDict
sed -i -E 's/^(\s*momentumAccommodationCoeff\s*)[0-9.]+;.*/\1 0;   \/\/ benchmark: recoil off/' system/controlDict
sed -i -E 's/^(\s*enable(PhaseExplosion|PlasmaModel|DropletBreakup)\s*)true;/\1false;/' system/controlDict
sed -i -E '/confinedVapourPocket/,/\}/ s/^(\s*enabled\s+)true;/\1false;   \/\/ benchmark: vapour pocket off/' system/controlDict
sed -i -E 's/^(\s*momentumBudget\s+)true;/\1false;/' system/fvSolution
sed -i -E "s/^endTime .*/endTime         $TEND;/; s/^maxDeltaT .*/maxDeltaT       $DTMAX;/; s/^writeInterval .*/writeInterval   $(python3 -c "print('%g'%($TEND/20))");/" system/controlDict
source /usr/lib/openfoam/openfoam2406/etc/bashrc >/dev/null 2>&1 || true
blockMesh > log.blockMesh 2>&1
rm -rf 0 && cp -r 0.orig 0
setExprFields > log.setExprFields 2>&1
echo "built $D: N=$N x $((2*N)) cells, dx=$DX nm, aspect=$ASP"
