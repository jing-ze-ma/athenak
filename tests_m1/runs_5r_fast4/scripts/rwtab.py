"""rwtab.py <dir with k*/ subdirs>: L1(rho), L1(E), Gamma error per case, nt, k; orders."""
import glob, json, os, sys, collections
R = sys.argv[1]
tab = collections.defaultdict(dict)
for f in glob.glob(R + '/k*/*/results.json'):
    d = json.load(open(f))
    c = d['case']
    k = os.path.basename(os.path.dirname(os.path.dirname(f)))
    ref = d['refs'][d['ref_key']]
    tab[(k, c['prat'], c['tau'])][c['nt']] = (d['l1_rho'], d['l1_E'],
                                               d['w_im'] - ref[1], d['w_re'] - ref[0])
for key in sorted(tab):
    row = tab[key]
    s = f'{key[0]:3s} P{key[1]:<5g} tau{key[2]:<6g}|'
    prev = None
    for nt in sorted(row):
        r, e, dg, dw = row[nt]
        o = '' if prev is None else f'({__import__("math").log2(prev/r):.2f})'
        s += f' {nt}: rho {r:.2e}{o} E {e:.2e} dG {dg:+.1e} |'
        prev = r
    print(s)
