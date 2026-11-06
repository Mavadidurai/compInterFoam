"""python3 fit_osc.py <case...>: fit (b-a)/R(t) = c + exp(-beta t) (A cos wt + B sin wt) by grid search over (w, beta)
with linear least squares for (c, A, B). Compares with Lamb's l = 2 result (low-viscosity limit)."""
import sys, numpy as np, math
sig=1.64; rho=4515; mu=2.35e-3; R=0.2e-6; nu=mu/rho
w0=math.sqrt(8*sig/(rho*R**3)); beta0=5*nu/R**2; wd=math.sqrt(w0**2-beta0**2)
print('Lamb l=2: omega0 = %.4e rad/s, damped omega = %.4e rad/s (period %.2f ns), beta = %.3e 1/s (tau %.1f ns)'%(w0,wd,2*math.pi/wd*1e9,beta0,1e9/beta0))
for c in sys.argv[1:]:
    D=np.atleast_2d(np.loadtxt(c+'/drop_metrics.csv',delimiter=',')); t=D[:,0]; y=(D[:,4]-D[:,5])/R
    best=None
    for w in np.linspace(0.5*wd,1.5*wd,801):
        for bt in np.linspace(0,4*beta0,161):
            M=np.c_[np.ones_like(t),np.exp(-bt*t)*np.cos(w*t),np.exp(-bt*t)*np.sin(w*t)]
            coef,res,_,_=np.linalg.lstsq(M,y,rcond=None); r=((M@coef-y)**2).sum()
            if best is None or r<best[0]: best=(r,w,bt,coef)
    r,w,bt,coef=best; amp=math.hypot(coef[1],coef[2])
    print('%-12s omega = %.4e rad/s (%+.1f %% vs Lamb), period %.2f ns; beta = %.3e 1/s (%+.0f %%); amplitude %.3f; rms residual %.4f; %d samples'%(c,w,100*(w/wd-1),2*math.pi/w*1e9,bt,100*(bt/beta0-1),amp,math.sqrt(r/len(t)),len(t)))
