"""Which force changes the momentum of the donor metal: force fields integrated with three weightings.

    xvfb-run -a pvpython tools/force_decomposition.py <case> <t_min [s]> <t_max [s]> <stride> <out.csv>

Needs a run with momentumBudget true, which writes budgetPressureForce (-grad p_rgh, reconstructed
from the face flux of the pressure solve) and budgetSurfaceTensionForce (the reconstructed CSF
face force) at every write. The pocket and recoil force density is taken as recoilPressure * dalpha/dy.

For every written time it integrates the y-component over the metal above the receiver (y > 8.05 um),
full circle, in units of 1e-13 kg m/s per ns, weighted
  _a  by alpha (as the logged momentumBudget lines), and
  _m  by the metal mass fraction alpha*rho_m/rho (the share of the cell force that acts on metal).
Columns Fm_*: total of the three forces, mass-weighted, by alpha class. ST_cap / ST_rest_*: the
mass-weighted surface tension on the cap (r < 1.5 um, y < 19.5 um) and on the rest, split at alpha 0.5.
ST_unweighted_all: the surface-tension force summed over all cells. For any set of closed interfaces
it is zero physically; a non-zero value is a numerical net force of the surface-tension discretisation.

Integrating a column over time and comparing with the change of the film momentum
(tools/mesh_study_fields.py) shows which force carries the change. The forces flicker from step to
step, so the time integral of written snapshots is uncertain by roughly +-0.1e-13 kg m/s per 3.6 ns
(docs/INVESTIGATION_LOG.md section 14).
"""
import os, sys, numpy as np
from paraview.simple import OpenFOAMReader, MergeBlocks, CellCenters, CellSize, Gradient, servermanager
from vtk.util.numpy_support import vtk_to_numpy as v2n
case=sys.argv[1]; T0,T1,STRIDE=float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]); f=os.path.join(case,'case.foam'); open(f,'a').close()
r=OpenFOAMReader(FileName=f); r.MeshRegions=['internalMesh']; r.CaseType='Decomposed Case'
r.CellArrays=['alpha.metal','U','p_rgh','recoilPressure','budgetPressureForce','budgetSurfaceTensionForce']; r.UpdatePipelineInformation()
mb=MergeBlocks(Input=r)
gr=Gradient(Input=mb); gr.ScalarArray=['CELLS','alpha.metal']; gr.ResultArrayName='ga'
cc=CellCenters(Input=CellSize(Input=gr)); W=72.0; RHO=4506.0; RG=1.1614
ts=[t for t in r.TimestepValues if T0<t<T1][::STRIDE]
rows=[]
for t in ts:
    cc.UpdatePipeline(t); d=servermanager.Fetch(cc); P=d.GetPointData(); g=lambda n: v2n(P.GetArray(n))
    X=v2n(d.GetPoints().GetData()); a=np.clip(g('alpha.metal'),0,1); V=g('Volume'); y=X[:,1]; rr=np.hypot(X[:,0],X[:,2])
    fp=g('budgetPressureForce')[:,1]; fs=g('budgetSurfaceTensionForce')[:,1]; pr=g('recoilPressure'); ga=g('ga')[:,1]
    frc=pr*ga   # recoil/pocket CSF force density, y
    m=y>8.05e-6
    wa=a; wm=a*RHO/(a*RHO+(1-a)*RG); wi=(a>1e-3).astype(float)
    def S_(fd,w,k=m): return (w*fd*V)[k].sum()*W*1e13*1e-9   # 1e-13 kg m/s per ns
    row=[t*1e9]
    for fd in (fp,fs,frc):
        row += [S_(fd,wa),S_(fd,wm)]
    # where is the mass-weighted pressure force: by alpha class and by height
    for lo,hi in ((0.5,1.01),(0.01,0.5),(0,0.01)):
        k=m&(a>=lo)&(a<hi); row.append(S_(fp+fs+frc,wm,k))
    row.append(S_(frc,wm,m&(y<19.5e-6)))   # pocket force on metal that has left the donor level
    cap=m&(rr<1.5e-6)&(y<19.5e-6); rest=m&~cap
    row += [S_(fs,wm,cap), S_(fs,wm,rest&(a<0.5)), S_(fs,wm,rest&(a>=0.5)), (fs*V)[m].sum()*W*1e13*1e-9, (fs*V).sum()*W*1e13*1e-9]
    rows.append(row)
    print(('t=%.3f  P: a %.3f m %.3f | ST: a %.3f m %.3f | pocket: a %.3f m %.3f | (P+ST+pk, mass-w) a>.5 %.3f  .01-.5 %.3f  <.01 %.3f | pocket on metal below 19.5um %.3f')%tuple(row[:11]) + '  ST cap %.3f rest-smeared %.3f rest-bulk %.3f unweighted film %.3f all %.3f'%tuple(row[11:]),flush=True)

np.savetxt(sys.argv[5],np.array(rows),delimiter=',',comments='',header='t_ns,P_a,P_m,ST_a,ST_m,pk_a,pk_m,Fm_a05,Fm_a001_05,Fm_a0_001,pk_below195,ST_cap,ST_rest_smeared,ST_rest_bulk,ST_unweighted_film,ST_unweighted_all')
