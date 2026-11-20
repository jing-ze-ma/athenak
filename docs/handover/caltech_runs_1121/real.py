# real.py RUNDIR... : per m1 dump in RUNDIR/a/bin, counts of |F_r|/(cE) > 1.05 / 1.2 / 1.5 in the top 120 radial shells (+ max)
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import bin_convert as bc
C=2.99792458e10; NT=120
for a in sys.argv[1:]:
    out=[]
    for f in sorted(glob.glob(os.path.join(a,'a','bin','*.m1.*.bin'))):
        d=bc.read_binary(f); md=d['mb_data']
        E=np.stack([np.asarray(md['m1_e'][m]) for m in range(d['n_mbs'])])[...,-NT:]
        F=np.stack([np.asarray(md['m1_f1'][m]) for m in range(d['n_mbs'])])[...,-NT:]
        r=np.abs(F/(C*E)); out.append('t%.0f %d/%d/%d max %.2f'%(d['time'],(r>1.05).sum(),(r>1.2).sum(),(r>1.5).sum(),r.max()))
    print(os.path.basename(a.rstrip('/')),' | '.join(out))
