import sys, glob, numpy as np
sys.path.insert(0,'/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'); import dhjcs
W='/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad/wbpair'
for d in sorted(glob.glob(W+'/*')):
    fs=sorted(glob.glob(d+'/bin/*.bin')); h=glob.glob(d+'/*.hst')
    print('==',d.split('/')[-1], 'last log:', open(d+'/out.txt').read().strip().split('\n')[-1][:60])
    if h:
        a=np.loadtxt(h[0]); print('   hst: t0 mass %.6e -> tN mass %.6e (%.2e); 1-KE: %.3e -> %.3e' % (a[0,2],a[-1,2],a[-1,2]/a[0,2]-1,a[0,7],a[-1,7]))
    for f in fs[-2:]:
        raw=dhjcs.bin_convert.read_binary(f); mb=raw['mb_data']; vr=np.asarray(mb['velx']); rho=np.asarray(mb['dens']); vt=np.sqrt(np.asarray(mb['vely'])**2+np.asarray(mb['velz'])**2)
        ni=vr.shape[-1]; q=lambda a,s: np.sqrt((a[...,s]**2).mean())
        print(f'   t={raw["time"]:.0f}s max|v_r|={np.abs(vr).max()/1e5:.3e} km/s rms v_r: bottom(i<8) {q(vr,slice(0,8))/1e5:.2e} mid {q(vr,slice(ni//2-4,ni//2+4))/1e5:.2e} top(i>=ni-8) {q(vr,slice(ni-8,ni))/1e5:.2e} | max|v_t|={vt.max()/1e5:.2e} | rho_min={rho.min():.2e}')
