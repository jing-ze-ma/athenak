#!/usr/bin/env python3
"""zzzzc5 item 4: mflux.py OUT DUMP ... : shell-mean 4 pi r^2 <rho v_r> (wedge solid-angle mean scaled to 4 pi)."""
import sys, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
out = sys.argv[1]
for f in sys.argv[2:]:
    h = bc.read_binary(f); n1 = h['mb_data']['dens'][0].shape[2]
    sr = np.zeros(n1); sm = np.zeros(n1); sv = np.zeros(n1); so = np.zeros(n1); rr = np.asarray(h['mb_x1v'][0], float); rf = np.asarray(h['mb_x1f'][0], float); om = 0; vmax = np.full(n1, -9e99)
    for b in range(h['n_mbs']):
        tf = np.asarray(h['mb_x2f'][b], float); pf = np.asarray(h['mb_x3f'][b], float)
        w = np.diff(pf)[:, None] * (np.cos(tf[:-1]) - np.cos(tf[1:]))[None, :]
        rho = np.asarray(h['mb_data']['dens'][b], float); v = np.asarray(h['mb_data']['velx'][b], float)
        W = w[:, :, None]; om += w.sum()
        sr += (W * rho * v).sum((0, 1)); so += (W * rho * np.maximum(v, 0)).sum((0, 1)); sm += (W * rho).sum((0, 1)); sv += (W * v * rho).sum((0, 1)); vmax = np.maximum(vmax, v.max((0, 1)))
    M = 4 * np.pi * rr ** 2 * sr / om; Mo = 4 * np.pi * rr ** 2 * so / om
    q = lambda x: int(np.argmax(x))
    i = q(M); o = q(Mo); j = q(-M); k = q(np.where(rr > 1.02 * rr[0], sv / sm, -9e99)); s1 = 'FILE %s t %.5e cyc %d' % (f.split('jingze/')[1], h['time'], h['cycle'])
    u = np.where(rr > 4e12, M, -9e99); iu = int(np.argmax(u))
    s1 += ' | r > 4e12: max net 4 pi r^2 <rho v_r> = %.4e g/s at r %.4e cm (<rho> %.3e, <v_r>_rho %.3e)' % (M[iu], rr[iu], sm[iu] / om, sv[iu] / sm[iu])
    open(out, 'a').write(s1 + ' | max net 4 pi r^2 <rho v_r> = %.4e g/s at r %.4e cm (shell <rho> %.3e, <v_r>_rho %.3e cm/s); max of outgoing part %.4e at r %.4e (cell v_r max there %.3e); most negative net %.4e at r %.4e; shell-mean max <v_r>_rho %.3e cm/s at r %.4e; net at r_top %.3e; ' % (M[i], rr[i], sm[i] / om, sv[i] / sm[i], Mo[o], rr[o], vmax[o], M[j], rr[j], sv[k] / sm[k], rr[k], M[-1]) + '\n')
    print('done', f, flush=True)
