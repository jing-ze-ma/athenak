# rstcmp.py <runs/JOBID>: restart checks for berx / rs10 / rs0 (TASK-2026-10-11-deltaai-agcarA-icrad-restart)
import sys, os, glob, hashlib, re, warnings
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, '/u/jma20/ATHENAK/builds/wt/hes_gpu/vis/python')
import bin_convert as bc
S = sys.argv[1]; B = f'{S}/berx'

def hrows(d, f='agcar3d.user.hst'):
    rows = {}
    for l in open(f'{d}/{f}'):
        if l.startswith('#'): continue
        rows.setdefault(l.split()[0], l)   # key: time string (exact text), first occurrence
    return rows

def tail_md5(f):   # md5 of the rst after the input-parameter text (which carries the command-line nlim)
    b = open(f, 'rb').read(); i = b.find(b'<par_end>'); return hashlib.md5(b[i:]).hexdigest(), len(b) - i

def bins_at(d, cyc):
    out = {}
    for f in glob.glob(f'{d}/bin/*.bin'):
        h = bc.read_binary(f)
        if h['cycle'] == cyc: out.setdefault(os.path.basename(f).split('.')[1], h)
    return out

for rs, c0, c1, rstb in (('rs10', 10, 20, 'agcar3d.00002.rst'), ('rs0', 0, 10, 'agcar3d.00001.rst')):
    d = f'{S}/{rs}'
    if not os.path.exists(f'{d}/out.txt'): print(f'{rs}: not run'); continue
    o = open(f'{d}/out.txt', errors='replace').read()
    icl = [l for l in o.splitlines() if 'ic_rad' in l]
    print(f'== {rs} (restart from cycle {c0}, to {c1})')
    print('  ic_rad lines in the restart log:'); [print('    ' + l) for l in icl]
    rerun = any(re.search(r'ic_rad (relax +\d|fs\d|fs_end|fs:|done)', l) for l in icl)
    nc = re.findall(r'NON-CONVERGED=(\S+)', o)[-1:]
    print(f'  IC step re-run on restart: {"YES" if rerun else "no"}; NON-CONVERGED: {nc}')
    for f in ('agcar3d.user.hst', 'agcar3d.hydro.hst'):
        hb, hr = hrows(B, f), hrows(d, f)
        common = [t for t in hr if t in hb]
        same = [hb[t] == hr[t] for t in common]
        cyc = sorted(round(float(t) / 1239.996) for t in common)
        print(f'  {f}: {len(common)} rows matched by time (cycles {cyc[0] if cyc else "-"}..{cyc[-1] if cyc else "-"}), '
              f'{sum(same)} identical text; rows only in restart: {len(hr) - len(common)}')
    bb, br = bins_at(B, c1), bins_at(d, c1)
    for k in sorted(set(bb) & set(br)):
        eq = all(np.array_equal(np.asarray(bb[k]['mb_data'][v]), np.asarray(br[k]['mb_data'][v])) for v in bb[k]['var_names'])
        print(f'  bin {k} at cycle {c1}: arrays bitwise {eq}')
    if not (set(bb) & set(br)): print(f'  no bin at cycle {c1} in both (berx bins at cycles 0, 20)')
    rr = sorted(glob.glob(f'{d}/rst/*.rst'))
    for f in rr:
        m, n = tail_md5(f); mb, nb = tail_md5(f'{B}/rst/{rstb}')
        print(f'  rst {os.path.basename(f)} vs berx {rstb} (cycle {c1}): data after <par_end> identical {m == mb} ({n} / {nb} bytes)')
