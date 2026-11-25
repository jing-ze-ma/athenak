import glob, re, numpy as np, sys
sys.path.insert(0,'/resnick/groups/carnegie_poc/jingze/pages_c1_1124/src')
RP=9.554e11; Om=0.121229
z=np.load('ring.npz'); t=z['t']; r=z['r']; R=z['R']; ir=lambda x:int(np.argmin(abs(r-x))); i=ir(0.74); ln=np.log(R[:,i])
for a,b in ((3e4,1.1e5),(1.1e5,1.85e5),(1.95e5,2.6e5),(2.6e5,3.2e5),(3.2e5,4e5),(4e5,7.8e5)):
    m=(t>=a)&(t<=b); p=np.polyfit(t[m],ln[m],1); print('ring rms r0.74 %.2e-%.2e n %d e-fold %.3e s'%(a,b,m.sum(),1/p[0]), 'start %.2e end %.2e km/s'%(np.exp(ln[m][0])/1e5,np.exp(ln[m][-1])/1e5))
mr=z['mean']; print('shell mean v_r at 0.74 (t>6e5): max |mean| %.3f km/s; ratio to ring rms at end %.1e'%(abs(mr[t>6e5][:,i]).max()/1e5,abs(mr[-1,i])/R[-1,i]))
# energy of the deep layer against hst
nm=None; H=[]
for l in open('/resnick/groups/carnegie_poc/jingze/prod_c1_1124/m35mms/m35mms.user.hst'):
    if l.startswith('#'):
        if '[1]=' in l: nm=re.findall(r'\]=(\S+)',l)
        continue
    p=l.split()
    if len(p)==len(nm): H.append([float(x) for x in p])
H=np.array(H); g={n:H[:,k] for k,n in enumerate(nm)}
D=[]
for f in sorted(glob.glob('out/m_*.npz')):
    d=np.load(f); rr=d['r']; m=(rr<0.80*RP); Er=0.5*Om*np.trapezoid(rr[m]**2*d['rw'][m],rr[m]); D.append((float(d['t']),Er))
D=np.array(D); 
print('t, 0.5*int rho v_r^2 dV over r<0.80 R_ph (incl. shell mean) [erg], KEr_int (hst) [erg], share')
for tt in (3.6e5,4.5e5,5.4e5,6.0e5,6.6e5,7.2e5,7.7e5):
    j=np.argmin(abs(D[:,0]-tt)); k=np.argmin(abs(g['time']-D[j,0])); print('%.3e %.3e %.3e %.2f  KE_int %.3e'%(D[j,0],D[j,1],g['KEr_int'][k],D[j,1]/g['KEr_int'][k],g['KE_int'][k]))
k0=np.argmin(abs(g['time']-6.0e5)); print('KEr_int gain from 6.0e5 to 7.7e5: %.2e; deep layer gain: %.2e'%(g['KEr_int'][np.argmin(abs(g['time']-7.7e5))]-g['KEr_int'][k0], D[np.argmin(abs(D[:,0]-7.7e5)),1]-D[np.argmin(abs(D[:,0]-6.0e5)),1]))
# numbers
P=np.load('prof.npz'); x=P['x']; kh=22*2*np.pi/(0.34906585*P['r']); 
print('lateral wavelength at 0.74: %.3e cm = %.2f cells (dtheta %.4f rad); kh %.2e /cm; kh*Hp %.1f'%(2*np.pi/kh[np.argmin(abs(x-0.74))],256/22,0.34906585/256,kh[np.argmin(abs(x-0.74))],kh[np.argmin(abs(x-0.74))]*P['Hp'][np.argmin(abs(x-0.74))]))
for xx in (0.72,0.74,0.76,0.78):
    j=np.argmin(abs(x-xx)); print('r %.2f N %.2e 1/(chi kh^2) %.2e s  Hp^2/chi %.2e s  Hp/cs %.2e s  2pi/N %.2e s  cs*kh %.2e rad/s (acoustic) acoustic cutoff cs/2Hp %.2e  ; kT %.2f kr %.2f  kT+kr/0.39 %.2f'%(x[j],np.sqrt(P['N2'][j]),1/(P['chi'][j]*kh[j]**2),P['Hp'][j]**2/P['chi'][j],P['Hp'][j]/P['cs'][j],2*np.pi/np.sqrt(P['N2'][j]),P['cs'][j]*kh[j],P['cs'][j]/(2*P['Hp'][j]),P['kT'][j],P['kr'][j],P['kT'][j]+P['kr'][j]/0.39))
