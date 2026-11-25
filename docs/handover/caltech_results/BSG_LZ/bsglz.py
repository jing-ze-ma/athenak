#!/usr/bin/env python3
"""zzzzd5 items 1-3 (READ ONLY).  user.hst Lz_tot / Lz_hflx per 500 cycles + MLATDIAG lines."""
import re, sys, glob, os, numpy as np
J = '/resnick/groups/carnegie_poc/jingze/'
def rd(f):
    names = None; rows = []
    for ln in open(f, errors='replace'):
        if ln.startswith('#'):
            if '[1]=' in ln: names = [m.split('=', 1)[1] for m in re.findall(r'\[\d+\]=\S+', ln)]
        elif names and len(ln.split()) == len(names): rows.append([float(x) for x in ln.split()])
    a = np.array(rows); t = a[:, 0]; keep = t < np.minimum.accumulate(np.append(t[1:], np.inf)[::-1])[::-1]
    return names, a[keep], len(a)
for lab, d, logs in (('c1 (95813cf85)', 'prod_c1_1124/bsg', None), ('T2 rc4v (262affd6f)', 'prod_rc4v_1124/bsg', None), ('rc3 (af836f5a3, chunk 4 dflfix 5b6159ad0)', 'prod_rc3_1124/bsg', None)):
    n, a, raw = rd(J + d + '/bsg3d.user.hst'); t = a[:, 0]; lz = a[:, n.index('Lz_tot')]; lf = a[:, n.index('Lz_hflx')]
    # cycle of every row: from the log lines cycle=.. time=..; rows are per cycle in the fresh c1 and T2 runs
    ct = []
    for f in sorted(glob.glob(J + d + '/log.[0-9]*'), key=os.path.getmtime):
        for ln in open(f, errors='replace'):
            if 'cycle=' in ln and 'time=' in ln and 'elapsed' in ln:
                m = re.search(r'cycle=(\d+) time=([0-9.eE+]+)', ln)
                if m: ct.append((float(m.group(2)), int(m.group(1))))
    ct = np.array(sorted(set(ct))); cyc = np.interp(t, ct[:, 0], ct[:, 1])
    print('\n## %s: rows %d (raw %d), t %.4e-%.4e, log cycle range %d-%d; Lz_tot first %.4e last %.4e, max |Lz_tot| %.4e' % (lab, len(a), raw, t[0], t[-1], ct[0, 1], ct[-1, 1], lz[0], lz[-1], abs(lz).max()))
    print('  cycles          Lz_tot start   end           change        change/500 cyc  Lz_hflx mean   residual = dLz + integral(Lz_hflx dt)   [g cm2/s; "Lz_hflx" as the hst column]')
    edges = np.arange(0, cyc[-1] + 500, 500)
    for lo in edges[:-1]:
        s = np.where((cyc >= lo) & (cyc < lo + 500))[0]
        if len(s) < 5: continue
        i, j = s[0], s[-1]; dt = np.diff(t[s]); fx = 0.5 * (lf[s][1:] + lf[s][:-1])
        print('  %6d-%-6d  %13.5e %13.5e %12.4e  %12.4e  %12.4e   %12.4e' % (lo, lo + 500, lz[i], lz[j], lz[j] - lz[i], (lz[j] - lz[i]) / max(cyc[j] - cyc[i], 1) * 500, lf[s].mean(), (lz[j] - lz[i]) + (fx * dt).sum()))
    for lo, hi in ((0, 5000), (5000, 6000), (6000, 7000), (7000, 8000), (8000, 9000), (9000, 10000)):
        s = np.where((cyc >= lo) & (cyc < hi))[0]
        if len(s) < 20: continue
        k = (cyc[s[-1]] - cyc[s[0]]) / 500
        print('  window %d-%d: Lz change per 500 cycles mean %.4e, |change| mean of 500-cycle windows %s' % (lo, hi, (lz[s[-1]] - lz[s[0]]) / max(k, 1e-9), ''))
    L = []
    for f in sorted(glob.glob(J + d + '/log.[0-9]*'), key=os.path.getmtime):
        for ln in open(f, errors='replace'):
            if ln.startswith('MLATDIAG'):
                w = ln.split(); m = re.search(r'MLATDIAG (\S+) cyc (\d+) t (\S+) dt (\S+) M (\S+) E (\S+) P (\S+) (\S+) (\S+) L (\S+) (\S+) (\S+) nfb (\d+) chk (\S+) (\S+) (\S+) (\S+)', ln)
                if m: L.append((m.group(1), int(m.group(2))) + tuple(float(x) for x in m.groups()[2:]))
    print('  MLATDIAG lines %d (tag cyc t L_z nfb chk_M chk_E chk_P chk_L); every line:' % len(L))
    for r in L: print('    %s c%-6d t %.5e Lx %.5e Ly %.5e Lz %.5e nfb %d chk %.1e %.1e %.1e %.1e' % (r[0], r[1], r[2], r[9], r[10], r[11], r[12], r[13], r[14], r[15], r[16]))
