from lib import *
from scipy.optimize import curve_fit
out=open('results4.txt','w')
def pr(*a):
    print(*a,flush=True); print(*a,file=out,flush=True)
PER=53.; TE=573.; ts=np.load('timescales.npz'); tth=ts['tth']; Rph=float(ts['Rph'])
z=np.load('hst_usr.npz'); n=list(z['names']); a=z['a']; t=a[:,0]/86400; Lin=a[-1,n.index('L_in')]
W7=[(TE-(k+1)*PER,TE-k*PER) for k in range(7)][::-1]; xc=np.array([.5*(u+v) for u,v in W7]); xs=xc*86400
wmn=lambda q,T=tp:np.array([q[(T>=u)&(T<v)].mean(0) for u,v in W7])
# 1. continuity vs dump mean mass flux
M=np.cumsum(Q['rho']*dV,axis=1)
Mw=wmn(M); Pc4=-np.polyfit(xs[-4:],Mw[-4:],1)[0]; Pd4=wmn(Q['Pm'])[-4:].mean(0)
RR=[9,10,12,15,17,20,25,30,40,50,60,70,80,90,95]; ii=[int(np.argmin(abs(r/RS-x))) for x in RR]
pr('MASS FLUX CHECK (1e22 g/s, 361-573): r, Pm(cell-centred rho v from dumps), -dM(<r)/dt (continuity)')
for R,i in zip(RR,ii): pr(' %5.1f %+7.2f %+7.2f'%(R,Pd4[i]/1e22,Pc4[i]/1e22))
# 2. Eulerian: L_in - dE(<r)/dt vs turbulent + radiative (cell-centred), with mean-flow term from continuity
Ecum=np.cumsum(Q['Etotd']*dV,axis=1); Ew=wmn(Ecum)
dE4=np.polyfit(xs[-4:],Ew[-4:],1)[0]; 
H=(5/3*Q['eint']+4/3*Q['E']+Q['ke'])/Q['rho']+Phi
Hw=wmn(H)[-4:].mean(0)
Lst=wmn(Q['Lrcom'])+wmn(Q['Lcv'])
Lnet=Lin-dE4; Lmf=Pc4*Hw
Lcheck=Lnet-Lmf
pr('\nFLUX CLOSURE (1e38, 361-573). A = L_in - dE(<r)/dt - Pc*H: conservative-budget flux in the turbulent+radiative channel; B = Lrcom+Lcv from dumps. (E incl KE+potential)')
pr(' r/Rs  A     B    B-A   Lrcom Lcv  B_std(4w)  Fconv/L_B  drift_B[%/100d,4w] SE')
for R,i in zip(RR,ii):
    s=fit(Lst[-4:,i],xc[-4:])
    pr(' %5.1f %6.3f %6.3f %+6.3f %6.3f %6.3f %6.3f %6.3f %+6.2f %5.2f'%(R,Lcheck[i]/1e38,Lst[-4:,i].mean()/1e38,(Lst[-4:,i].mean()-Lcheck[i])/1e38,wmn(Q['Lrcom'])[-4:,i].mean()/1e38,wmn(Q['Lcv'])[-4:,i].mean()/1e38,Lst[-4:,i].std()/1e38,wmn(Q['Lcv'])[-4:,i].mean()/Lst[-4:,i].mean(),1e4*s[0]/s[3],1e4*s[1]/s[3]))
c=(r>25*RS)&(r<90*RS); pr('B over 25-90 Rs: mean %.3f, std over radius %.3f (1e38); L_in %.3f, L_top(361-573) %.3f'%(Lst[-4:,c].mean()/1e38,Lst[-4:,c].mean(0).std()/1e38,Lin/1e38,a[(t>=361)&(t<573),n.index('L_top')].mean()/1e38))
# per-window Lrcom,Lcv at mid-CZ
for k in ('Lrcom','Lcv','Lconv','Lkin'):
    q=wmn(Q[k]); pr(k,'/L windows @25,40,60,80 Rs:',[ [round(q[w,int(np.argmin(abs(r/RS-x)))]/1.38e38,3) for x in (25,40,60,80)] for w in range(3,7)])
# 3. Lagrangian: heat storage L(m)=L_in - dE(<m)/dt - P dV/dt
Mref=M[tp>361][:,r<=Rph][:,-1].mean()
Fr=np.array([.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95])
Pt=Q['E']/3+2/3*Q['eint']
rfo=rf[1:]
RM=np.zeros((len(tp),len(Fr)));EM=RM.copy();PM=RM.copy();UM=RM.copy()
Ucum=np.cumsum(Q['U']*dV,axis=1)
for j in range(len(tp)):
    RM[j]=np.interp(Fr*Mref,M[j],rfo); EM[j]=np.interp(RM[j],rfo,Ecum[j]); UM[j]=np.interp(RM[j],rfo,Ucum[j]); PM[j]=np.interp(RM[j],r,Pt[j])
Vw=wmn(4*np.pi/3*RM**3); Ew2=wmn(EM); Pw=wmn(PM)
pr('\nLAGRANGIAN HEAT STORAGE: Lq(m) = L_in - dE(<m)/dt - P dV_m/dt [1e38] (<<KH storage: L_in - Lq = int T ds/dt dm>>); slope on 4w (361-573) and 7w (202-573); SE; last two windows halves')
pr(' m  r/Rs  Lq4 SE4 | Lq7 SE7 | Storage4/L[%] Storage7/L[%]')
Stor=[]
for k,f in enumerate(Fr):
    o=[]
    for nw in (4,7):
        a1=fit(Ew2[-nw:,k],xs[-nw:]); b1=fit(Vw[-nw:,k],xs[-nw:]); pm=Pw[-nw:,k].mean()
        Lq=Lin-a1[0]-pm*b1[0]; se=np.hypot(a1[1],pm*b1[1]); o.append((Lq,se))
    Stor.append((f,RM[tp>361,k].mean()/RS,o[0][0],o[0][1],o[1][0],o[1][1]))
    pr('%.2f %6.1f %6.3f %5.3f | %6.3f %5.3f | %+6.2f %+6.2f'%(f,RM[tp>361,k].mean()/RS,o[0][0]/1e38,o[0][1]/1e38,o[1][0]/1e38,o[1][1]/1e38,100*(Lin-o[0][0])/Lin,100*(Lin-o[1][0])/Lin))
np.savez('lagr.npz',Stor=np.array(Stor),Lin=Lin,r=r,Lcheck=Lcheck,Lst=Lst,xc=xc,Lrcom=wmn(Q['Lrcom']),Lcv=wmn(Q['Lcv']),Fr=Fr,RM=RM,tp=tp)
# 4. exponential decay of the secular part on P-window means
def ex(x,qi,A,tau): return qi+A*np.exp(-(x-xc[0])/tau)
pr('\nEXP FIT q(t)=q_inf+A exp(-(t-228)/tau) on 7 window means (ln q, in %). tau[d] with fit SE; rem = offset of last window from q_inf [%]; noise = residual std [%] ')
for lab,k in (('rho','rho'),('T','Tgas')):
    pr(lab,' r tau SE A% rem% noise%')
    Y=wmn(Q[k])
    for R,i in zip(RR,ii):
        y=100*np.log(Y[:,i]); y=y-y[-1]
        try:
            p,cv=curve_fit(ex,xc,y,p0=(0,y[0],100),bounds=([-50,-100,10],[50,100,5000]),maxfev=5000); e=y-ex(xc,*p)
            pr(' %5.1f %6.0f %6.0f %+7.2f %+6.2f %5.2f  tth %.0f'%(R,p[2],np.sqrt(cv[2,2]),p[1],y[-1]-p[0],e.std(ddof=3),tth[i]))
        except Exception as ee: pr(' %5.1f fit failed %s'%(R,ee))
Y=wmn(RM); pr('r(m) tau')
for k,f in enumerate(Fr):
    y=100*np.log(Y[:,k]); y=y-y[-1]
    p,cv=curve_fit(ex,xc,y,p0=(0,y[0],100),bounds=([-50,-100,10],[50,100,5000]),maxfev=5000); e=y-ex(xc,*p)
    pr(' m %.2f r %.1f tau %.0f SE %.0f A %+.2f rem %+.2f noise %.2f'%(f,RM[tp>361,k].mean()/RS,p[2],np.sqrt(cv[2,2]),p[1],y[-1]-p[0],e.std(ddof=3)))
# 5. vs IC
pr('\nCOMPARE TO IC (t<2 d mean) : % change of last window (520-573) at fixed r')
ic=tp<2
for R,i in zip(RR,ii): pr(' %5.1f rho %+7.2f T %+7.2f'%(R,100*(wmn(Q['rho'])[-1,i]/Q['rho'][ic,i].mean()-1),100*(wmn(Q['Tgas'])[-1,i]/Q['Tgas'][ic,i].mean()-1)))
# 6. pulsation coherence in depth (sinusoid at P=53, 361-573)
pr('\nPULSATION in r(m): amplitude[%], phase[deg] for P=53 d on 361-573')
s=tp>361
for k,f in enumerate(Fr):
    A_=np.c_[np.cos(2*np.pi*tp[s]/53),np.sin(2*np.pi*tp[s]/53),np.ones(s.sum()),tp[s]-470]; co,*_=np.linalg.lstsq(A_,RM[s,k],rcond=None)
    pr(' m %.2f r %.1f amp %.3f%% phase %.0f'%(f,co[2]/RS,100*np.hypot(co[0],co[1])/co[2],np.degrees(np.arctan2(co[1],co[0]))))
# 7. wall rate
Hd=H_=np.load('dump_hdr.npy'); ti=Hd[:,1]/86400; mt=Hd[:,2]/86400
for span in (3,10,30):
    s=ti>ti[-1]-span
    # fit over recent span; rate in star-d per wall-d
    pr('wall rate last %d star-d: dump count %d, wall span %.2f d, %.1f star-d/wall-d'%(span,s.sum(),mt[s][-1]-mt[s][0],(ti[s][-1]-ti[s][0])/(mt[s][-1]-mt[s][0])))
