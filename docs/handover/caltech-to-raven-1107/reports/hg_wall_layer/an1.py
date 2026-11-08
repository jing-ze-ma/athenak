from lib import *
import sys
out=open('results1.txt','w')
def pr(*a):
    print(*a,flush=True); print(*a,file=out,flush=True)
z=np.load('hst_usr.npz'); n=list(z['names']); a=z['a']; t=a[:,0]/86400; gh=lambda k:a[:,n.index(k)]
F=np.load('hst_usr_full.npy'); tf=F[:,0]/86400; gf=lambda k:F[:,n.index(k)]
pr('stale reds dropped',stale.sum(),'of',len(stale),'range t',tp[0],tp[-1])
# event: strip
d=np.diff(gf('M_tot')); j=np.argsort(d)[:3]; pr('M_tot steps',[(round(tf[x],3),'%.3e'%d[x]) for x in j]); ib=j[0]; pr('Etot across step %.4e -> %.4e'%(gf('Etot')[ib],gf('Etot')[ib+1]))
t_s=tf[ib]
pr('red rho at r: t before/after strip %.2f'%t_s)
b=np.argmin(abs(tp-(t_s-0.4))); c=np.argmin(abs(tp-(t_s+0.4))); pr('dumps',tp[b],tp[c])
dm=(Q['rho'][b]-Q['rho'][c])*dV; cum=np.cumsum(dm[::-1])[::-1]
for R in (60,80,90,100,110,120,130,140,150,170,200,300):
    i=np.argmin(abs(r/RS-R)); pr(' r %4d Rsun rho %.3e -> %.3e  mass removed above r %.3e'%(R,Q['rho'][b][i],Q['rho'][c][i],cum[i]))
pr('Rph before/after',Q['Rph_med'][b]/RS,Q['Rph_med'][c]/RS)
s=(tf>526.5)&(tf<527.1); L=gf('L_top')[s]; i=np.argsort(L)[-3:]; pr('spike 526.8 full-res: t,L_top,L_mid',[(round(tf[s][k],4),round(L[k]/1e38,3),round(gf('L_mid')[s][k]/1e38,3)) for k in i],'n rows above 1.42e38',(L>1.42e38).sum(),'span d',np.ptp(tf[s][L>1.42e38]))
# period
s=t>325; dt=np.median(np.diff(t[s])); ti=np.arange(t[s][0],573,dt)
for nm,tt,yy in (('Rph',tp,Q['Rph_med']),('L_top',t,gh('L_top')),('Pm@60Rs',tp,Q['Pm'][:,np.argmin(abs(r/RS-60))])):
    yi=np.interp(ti,tt,yy); yi=yi-np.polyval(np.polyfit(ti,yi,2),ti)
    f=np.linspace(1/120,1/5,2000); A=[abs(np.sum(yi*np.exp(-2j*np.pi*fq*ti))) for fq in f]; A=np.array(A)*2/len(ti); pk=[k for k in range(1,len(A)-1) if A[k]>A[k-1] and A[k]>A[k+1]]; pk=sorted(pk,key=lambda k:-A[k])[:5]
    pr('period',nm,[(round(1/f[k],2),'%.3e'%A[k]) for k in pk])
