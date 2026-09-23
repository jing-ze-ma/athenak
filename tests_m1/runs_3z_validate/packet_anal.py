#!/usr/bin/env python3
"""Free-streaming Gaussian packet (problem/beam_packet): python3 packet_anal.py <run dirs>
Per run, final bin dump vs the first: energy change, peak E ratio, centroid travel / (c t),
RMS width PERPENDICULAR (transverse, the beam-width measure) and PARALLEL to the motion,
both relative to the initial (exact free streaming keeps all of them), and the negative-E /
|F| > cE cell counts.  Periodic box: coordinates are unwrapped about the expected centre."""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/tests_m1')
import common  # noqa: E402

N1 = N2 = 1.0 / np.sqrt(2.0)
X0 = Y0 = 0.3
C = 1.0


def mom(d, xc, yc):
    e = np.asarray(d.var('m1_e'))[0] - 1.0e-15
    f1 = np.asarray(d.var('m1_f1'))[0]
    f2 = np.asarray(d.var('m1_f2'))[0]
    X, Y = np.meshgrid(d.x1, d.x2)
    dx = (X - xc + 0.5) % 1.0 - 0.5
    dy = (Y - yc + 0.5) % 1.0 - 0.5
    sp = dx * N1 + dy * N2
    sq = -dx * N2 + dy * N1
    w = e.sum()
    cp, cq = (e * sp).sum() / w, (e * sq).sum() / w
    vp = (e * (sp - cp) ** 2).sum() / w
    vq = (e * (sq - cq) ** 2).sum() / w
    fr = np.sqrt(f1 ** 2 + f2 ** 2) / np.maximum(C * (e + 1e-15), 1e-300)
    return dict(w=w, peak=e.max(), cp=cp, cq=cq, sp=np.sqrt(vp), sq=np.sqrt(vq),
                neg=int((e + 1e-15 <= 0).sum()), fmax=float(fr.max()),
                nf=int((fr > 1.0 + 1e-10).sum()))


for rd in sys.argv[1:]:
    fs = sorted(glob.glob(os.path.join(rd, 'bin', '*.bin')))
    if len(fs) < 2:
        print('%-28s no dumps' % os.path.basename(rd))
        continue
    d0, d1 = common.load_dump(fs[0]), common.load_dump(fs[-1])
    t = d1.time
    xc, yc = X0 + C * t * N1, Y0 + C * t * N2
    m0, m1 = mom(d0, X0, Y0), mom(d1, xc, yc)
    log = open(os.path.join(rd, 'run.log')).read()
    nc = re.findall(r'NON-CONVERGED=([0-9.eE+-]+)', log)
    fb = re.findall(r'stage fallbacks=([0-9.eE+-]+)', log)
    cyc = re.findall(r'cycle=(\d+)', log)
    print('%-26s t=%.3f nx=%d dt=%.3e | dE/E %+.2e peak %.3f | travel/(ct) %.4f '
          'dq %+.4f | sig_perp %.3f sig_par %.3f (x initial) | negE %d |f|>1 %d '
          'fmax %.4f | NC %s fallb %s'
          % (os.path.basename(rd), t, d1.x1.size, t / max(int(cyc[-1]), 1),
             m1['w'] / m0['w'] - 1.0, m1['peak'] / m0['peak'], 1.0 + m1['cp'] / (C * t),
             m1['cq'], m1['sq'] / m0['sq'], m1['sp'] / m0['sp'], m1['neg'], m1['nf'],
             m1['fmax'], nc[-1] if nc else '-', fb[-1] if fb else '-'))
