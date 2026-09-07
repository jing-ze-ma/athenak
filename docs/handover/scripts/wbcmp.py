import sys, glob, numpy as np
sys.path.insert(0,'/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad'); import dhjcs
W='/tmp/claude-28895/-viper-u2-jinma-ATHENAK-athenak/967431f9-6fb9-4e65-a716-2a3c3456d615/scratchpad/wbpair'
bands=[(0,16),(16,64),(64,112),(112,128)]
def stats(f):
    raw=dhjcs.bin_convert.read_binary(f); vr=np.asarray(raw['mb_data']['velx'])/1e5
    return raw['time'], [np.sqrt((vr[...,a:b]**2).mean()) for a,b in bands], np.abs(vr).max()
for pair in ['cs','sp','csmhd']:
    print(f'== {pair}: rms v_r [km/s] in radial bands i<16 | 16-64 | 64-112 | >=112, max|v_r|   (off / on)')
    fo={round(stats(f)[0]/1000): f for f in sorted(glob.glob(f'{W}/{pair}_off/bin/*.bin'))}
    fn={round(stats(f)[0]/1000): f for f in sorted(glob.glob(f'{W}/{pair}_on/bin/*.bin'))}
    for k in sorted(set(fo)&set(fn)):
        if k==0: continue
        t,so,mo=stats(fo[k]); _,sn,mn=stats(fn[k])
        print(f'   t={t:6.0f}s  off: '+' '.join(f'{v:.2e}' for v in so)+f'  max {mo:.2e}   |  on: '+' '.join(f'{v:.2e}' for v in sn)+f'  max {mn:.2e}   ratio on/off: '+' '.join(f'{b/a if a>0 else 0:.2f}' for a,b in zip(so,sn)))
    try:
        ho=np.loadtxt(glob.glob(f'{W}/{pair}_off/*.hst')[0]); hn=np.loadtxt(glob.glob(f'{W}/{pair}_on/*.hst')[0])
        n=min(len(ho),len(hn)); print('   hst 1-KE off/on at t=%.0f: %.3e / %.3e ; mass drift off/on: %.2e / %.2e'%(hn[n-1,0],ho[n-1,7],hn[n-1,7],ho[n-1,2]/ho[0,2]-1,hn[n-1,2]/hn[0,2]-1))
    except Exception as e: print('   hst:',e)
