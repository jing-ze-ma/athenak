import sys, glob, os, numpy as np
def last(d):
    f = sorted(glob.glob(d + "/tab/rw.m1.*.tab"))[-1]
    return np.loadtxt(f)[:, 3]   # m1_e
root = sys.argv[1]
groups = {}
for d in glob.glob(root + "/*"):
    b = os.path.basename(d)
    if not os.path.exists(d + "/results.json"): continue
    nt = int(b.split("_nt")[1].split("_")[0]); key = b.replace("_nt%d" % nt, "")
    groups.setdefault(key, {})[nt] = d
for key, g in sorted(groups.items()):
    nts = sorted(g); e = [last(g[n]) for n in nts]
    df = [np.abs(e[i] - e[i+1]).mean() for i in range(len(e)-1)]
    ords = [np.log2(df[i]/df[i+1]) for i in range(len(df)-1)]
    print(key.split("_amp")[0], "nts", nts, "diffs", ["%.3e" % x for x in df], "order", ["%.2f" % o for o in ords])
