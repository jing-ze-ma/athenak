import numpy as np,glob,os
edges=[0,100,200,400,600,800,850,900,1000,2000]
out=[]
for f in sorted(glob.glob('phase/*.imprint.npz')):
    z=np.load(f); r=z['r']; n=os.path.basename(f).replace('.bin.imprint.npz','')
    out.append('%s cycle %d t = %.3f d'%(n,z['cycle'],z['time']/86400))
    out.append('  r bin      max m4     r(max)  phi0[deg]  max m8     sym_rms    rms_all   sym/all*sqrt6  rms_panel min-max at r(max)')
    for lo,hi in zip(edges[:-1],edges[1:]):
        s=np.where((r>=lo)&(r<hi))[0]
        if len(s)==0: continue
        k=s[np.argmax(z['a4'][s])]
        ph=np.degrees(-np.angle(z['c4'][k])/4)%90
        kk=s[np.argmax(z['rms_all'][s])]
        out.append('  %4d-%-5d  %.3e %6.1f   %5.1f     %.3e  %.3e  %.3e  %.3f         %.2e-%.2e'%(lo,hi,z['a4'][k],r[k],ph,z['a8'][s].max(),z['rms_sym'][s].max(),z['rms_all'][s].max(),z['rms_sym'][kk]/z['rms_all'][kk]*6**.5,z['rms_panel'][:,k].min(),z['rms_panel'][:,k].max()))
open('phase_table.txt','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
b=np.load('phase/B_00008.bin.imprint.npz'); o=np.load('/resnick/groups/carnegie_poc/jingze/cshero_1124/B/run_cs/bin/rsg.hydro_w.00008.bin.imprint.npz')
print('B 00008 new vs old a4 max abs diff',abs(b['a4']-o['a4']).max())
m=lambda f:(lambda z:z['a4'][(z['r']>=200)&(z['r']<400)].max())(np.load(f))
for n in ('08362','08499','08636'): print('RSG %s / B 00008 m4(200-400) = %.1f'%(n,m('phase/RSG_%s.bin.imprint.npz'%n)/m('phase/B_00008.bin.imprint.npz')))
