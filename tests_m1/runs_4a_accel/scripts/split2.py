# usage: python3 split2.py <prof dir> <ncyc_total> <skip> <regex> : per-kernel ms/cycle for names matching regex
import csv, re, sys, collections, glob
sys.argv += []
d, ntot, skip, pat = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
exec(open('/viper/ptmp2/jinma/accel_0923/tools/split.py').read().split("CAT = [")[0].split("d, rk, ntot")[0])
def short(n):
    if 'M1PCRX<float>' in n: return 'M1PCRX<float>'
    if 'M1PCRX<double>' in n: return 'M1PCRX<double>'
    m = re.search(r'([A-Za-z_][\w:]*)\(\)[^:]*::\{lambda', n) or re.search(r'(\w[\w:]*)\([^()]*\)::\{lambda', n)
    s = m.group(1) if m else n[:60]
    lam = re.search(r'::\{lambda\([^)]*\)#(\d+)\}', n)
    return s + ('#' + lam.group(1) if lam else '')
kt = glob.glob(f'{d}/**/rank_0_kernel_trace.csv', recursive=True)[0]
rows = sorted((int(r['Start_Timestamp']), int(r['End_Timestamp']), short(r['Kernel_Name'])) for r in csv.DictReader(open(kt)))
T = [r[0] for r in rows if 'NewTimeStep' in r[2]]
per = max(1, round(len(T)/ntot)); t0, t1 = T[per*skip], T[-1]
w = [r for r in rows if t0 <= r[0] < t1]; ncyc = len([t for t in T if t0 <= t < t1])/per
byk = collections.Counter(); byn = collections.Counter()
for a, b, s in w:
    if re.search(pat, s): byk[s] += b - a; byn[s] += 1
f = 1e6*ncyc; tot = 0
for s, v in byk.most_common(60):
    tot += v; print(f'{v/f:7.3f} {byn[s]/ncyc:7.1f} {v/byn[s]/1e3:8.1f}  {s[:90]}')
print('total', tot/f)
