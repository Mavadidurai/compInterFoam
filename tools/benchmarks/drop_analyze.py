"""pvpython drop_analyze.py <case> [out.csv]
Per written time: metal volume (sum alpha*V), Laplace jump (mean p_rgh inside r<0.5R minus mean p_rgh in gas
beyond 2R), max|U| (spurious currents), drop half-extent along the axis (b) and at the equator (a) from the
alpha = 0.5 crossings.  R = 0.2 um drop centred at the origin."""
import sys, os, numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n
case=sys.argv[1]; out=sys.argv[2] if len(sys.argv)>2 else os.path.join(case,'drop_metrics.csv')
f=os.path.join(case,'case.foam'); open(f,'a').close()
r=OpenFOAMReader(FileName=f); r.MeshRegions=['internalMesh']; r.CaseType='Reconstructed Case'; r.SkipZeroTime=0
r.CellArrays=['alpha.metal','p_rgh','U']; r.UpdatePipelineInformation()
mb=MergeBlocks(Input=r); cs=CellSize(Input=mb); cs.ComputeVolume=1; cc=CellCenters(Input=cs)
R=0.2e-6; rows=[]
def crossing(s,a):
    o=np.argsort(s); s,a=s[o],a[o]; i=np.where((a[:-1]-0.5)*(a[1:]-0.5)<=0)[0]
    if len(i)==0: return np.nan
    i=i[-1]; return s[i]+(0.5-a[i])*(s[i+1]-s[i])/(a[i+1]-a[i])
for t in r.TimestepValues:
    cc.UpdatePipeline(t); d=servermanager.Fetch(cc); X=v2n(d.GetPoints().GetData()); P=d.GetPointData()
    a=v2n(P.GetArray('alpha.metal')); p=v2n(P.GetArray('p_rgh')); U=v2n(P.GetArray('U')); V=v2n(P.GetArray('Volume'))
    rad=np.hypot(X[:,0],X[:,1])
    inside=(rad<0.5*R)&(a>0.99); outside=(rad>2*R)&(a<0.01)
    dp=p[inside].mean()-p[outside].mean() if inside.any() and outside.any() else np.nan
    # axis column: smallest x; equator row: |y| smallest
    xmin=X[:,0].min(); ax=np.isclose(X[:,0],xmin,rtol=0,atol=1e-12)
    yabs=np.abs(X[:,1]); ymin=yabs.min(); eq=np.isclose(yabs,ymin,rtol=0,atol=1e-12)&(X[:,1]>0)
    b=crossing(X[ax,1][X[ax,1]>0],a[ax][X[ax,1]>0]); aeq=crossing(X[eq,0],a[eq])
    rows.append((t,(a*V).sum(),dp,np.linalg.norm(U,axis=1).max(),b,aeq,p[inside].std() if inside.any() else np.nan))
R_=np.array(rows); np.savetxt(out,R_,delimiter=',',header='t_s,metal_volume_m3,laplace_dp_Pa,max_U_ms,b_axis_m,a_equator_m,p_inside_std_Pa')
print('wrote',out,len(rows),'times')
