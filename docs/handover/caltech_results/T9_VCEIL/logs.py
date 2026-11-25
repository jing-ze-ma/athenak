#!/usr/bin/env python3
"""zzzzc4 item 4: counters from existing logs (plain read-only opens).  usage: logs.py DIR[:label] ...
t9l62.log event counters (per log interval, all c2p calls): eos_dfloor, eos_vceil, eos_tfloor, fofc, vceil_de.
floor_audit lines: cell-c2p, floored, skip(dfloor/vceil) are CUMULATIVE since the process start; increments are differenced per file
(a decrease = a new process).  skip(dfloor/vceil) mixes both floors (source: one counter DC_SKIP): NOT split."""
import os, re, sys, glob, numpy as np
WIN = (4.6e5, 5.7e5)
for a in sys.argv[1:]:
    d, _, lab = a.partition(':'); lab = lab or d
    print('\n## %s  (%s)' % (lab, d))
    # cycle -> time from the log lines
    ct = []
    for f in sorted(glob.glob(d + '/log.[0-9]*'), key=os.path.getmtime):
        L = open(f, errors='replace').read()
        ct += [(int(c), float(t)) for c, t in re.findall(r'cycle=(\d+) time=([0-9.eE+]+)', L)]
        fa = re.findall(r'floor_audit cycle (\d+) \[(update|refresh)\] cell-c2p ([0-9.eE+]+) .*?floored ([0-9.eE+]+) .*?skip\(dfloor/vceil\) ([0-9.eE+-]+) skip_de ([0-9.eE+-]+) skip_dm ([0-9.eE+-]+)', L)
        fc = re.findall(r'[^\n]*(?:floorcap|FloorGasCap)[^\n]*', L)
        if fa:
            by = {}
            for c, k, n, fl, sk, sde, sdm in fa:
                by.setdefault(k, []).append((int(c), float(n), float(fl), float(sk), float(sde)))
            for k, v in by.items():
                v = np.array(v); dn = np.diff(v[:, 1]); dfl = np.diff(v[:, 2]); ds = np.diff(v[:, 3]); cyc = v[1:, 0]; ok = (dn >= 0) & (ds >= 0)
                if ok.sum():
                    print('  %s [%s]: %d audit lines, cycles %d-%d; final cumulative cell-c2p %.4e floored %.4e skip %.4e (skip/cell-c2p %.3e, floored/cell-c2p %.3e, skip_de %.3e)' % (os.path.basename(f), k, len(v), v[0, 0], v[-1, 0], v[-1, 1], v[-1, 2], v[-1, 3], v[-1, 3] / max(v[-1, 1], 1), v[-1, 2] / max(v[-1, 1], 1), v[-1, 4]))
                    print('     per audit line: skip median %.3e p90 %.3e max %.3e at cycle %d; floored median %.3e max %.3e' % (np.median(ds[ok]), np.percentile(ds[ok], 90), ds[ok].max(), cyc[ok][np.argmax(ds[ok])], np.median(dfl[ok]), dfl[ok].max()))
        if fc:
            print('  %s: %d floorcap lines; last: %s' % (os.path.basename(f), len(fc), fc[-1].strip()[:330]))
    f = d + '/t9l62.log'
    if os.path.exists(f):
        hdr, rows = None, []
        for ln in open(f, errors='replace'):
            if ln.startswith('#  cycle'):
                hdr = ln[1:].split()
            elif not ln.startswith('#') and hdr and len(ln.split()) == len(hdr):
                rows.append([float(x) for x in ln.split()])
        if rows:
            R = np.array(rows); ix = {k: i for i, k in enumerate(hdr)}
            sel = {k: R[:, ix[k]] for k in ('eos_dfloor', 'eos_vceil', 'eos_tfloor', 'fofc', 'vceil_de', 'dfloor_ke') if k in ix}
            print('  t9l62.log: %d rows, cycles %d-%d' % (len(R), R[0, 0], R[-1, 0]))
            for k, v in sel.items():
                print('     %-11s sum %.4e  max per interval %.4e  intervals > 0: %d / %d' % (k, v.sum(), v.max(), (v > 0).sum(), len(v)))
            for c, vc, df in zip(R[:, 0], sel['eos_vceil'], sel['eos_dfloor']):
                if vc > 0:
                    print('     interval ending cycle %d: eos_vceil %d, eos_dfloor %d' % (c, vc, df))
    if ct:
        print('  cycle-time span in the logs: cycle %d t %.5e .. cycle %d t %.5e' % (ct[0][0], ct[0][1], ct[-1][0], ct[-1][1]))
