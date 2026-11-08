import numpy as np, os
RS=6.957e10; G=6.674e-8; GMC=3.6109195512e26; LIN=None
P_=np.load('prof.npz'); H=np.load('dump_hdr.npy')
RED='/resnick/groups/carnegie_poc/jingze/hgglob_page/out/data'
idx=P_['idx'].astype(int)
mt={int(h[0]):h[2] for h in H}
stale=np.array([os.path.getmtime(f'{RED}/red_{i:05d}.npz')<mt[i] for i in idx])
keep=~stale
tp=P_['t'][keep]/86400; r=P_['r']; rf=P_['rf'][:len(r)+1]
dV=4*np.pi/3*np.diff(rf**3)
Q={k:P_[k][keep] for k in 'rho mvr vrrms vtrms E Trad Tgas eint Lr Lconv Lconv_gas Lradadv Lrcom Lbulk Lradbulk Lkin dlrho_rms Rph_med Rph_p2 Rph_p98 Teff_med'.split()}
# gravity: IC column + core (he_gm_column), potential zero at the wall
m0=np.cumsum(P_['rho'][0]*dV)
g=(GMC+G*m0)/r**2
Phi=np.concatenate([[0],np.cumsum(0.5*(g[1:]+g[:-1])*np.diff(r))])
Q['Pm']=4*np.pi*r**2*Q['rho']*Q['mvr']                # mean mass flux [g/s]
Q['Lfl']=Q['Lconv_gas']+Q['Lkin']+Q['Lbulk']-Q['Lradbulk']+Q['Pm']*Phi     # gas energy flux incl. mean flow + potential
Q['LE']=Q['Lr']+Q['Lfl']                               # total conservative energy flux
Q['Lcv']=Q['Lconv']+Q['Lkin']                          # turbulent: gas+rad enthalpy + KE
Q['U']=Q['eint']+Q['E']
Q['ke']=0.5*Q['rho']*(Q['vrrms']**2+Q['mvr']**2+Q['vtrms']**2)
Q['Etotd']=Q['U']+Q['ke']+Q['rho']*Phi                 # energy density incl. potential
def wmean(q,a,b,T=None):
    T=tp if T is None else T
    s=(T>=a)&(T<b); return q[s].mean(0)
def fit(y,x=None):
    """OLS on independent window means: slope per unit x, SE, residual std"""
    n=len(y); x=np.arange(n,dtype=float) if x is None else x
    A=np.c_[x-x.mean(),np.ones(n)]; c,*_=np.linalg.lstsq(A,y,rcond=None); res=y-A@c
    s=np.sqrt((res**2).sum()/max(n-2,1)); se=s/np.sqrt(((x-x.mean())**2).sum())
    return c[0],se,s,c[1]
