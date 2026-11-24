#!/usr/bin/env python3
"""ROTALL ray-effect measure from one vet dump (+ optional hydro_w dump) of an sp wedge.
usage: raymeas.py <x.vet.NNNNN.bin> [<x.hydro_w.NNNNN.bin>] [fthin=0.5]
Per radial shell i: p(k) = sum_j F3 / sum_j |F1| over phi index k (ZAMSVPHI fprof form).
Thin shells: flux factor sum|F1| / sum E >= fthin x its max over shells. Prints the odd part
max_k |p(k) - p(K-1-k)|/2 (max and median over thin shells), the mean <p>, and, with the
hydro dump, v_phi rms over the thin shells and the mass-weighted v_phi rms of the whole mesh."""
import sys, types, numpy as np
sys.modules.setdefault('h5py', types.ModuleType('h5py'))
sys.path.insert(0, '/orion/ptmp/jinma/clean_1115/int2_tests/ana')   # Caltech: dir of bin_convert.py
import bin_convert as bc


def full(d, n):
    x2 = np.array(d['mb_x2f']); x3 = np.array(d['mb_x3f']); a = d['mb_data'][n]
    J = sorted(set(np.round(x2[:, 0], 12))); K = sorted(set(np.round(x3[:, 0], 12)))
    nk, nj, ni = np.array(a[0]).shape
    F = np.zeros((nk*len(K), nj*len(J), ni))
    for b in range(len(a)):
        j = J.index(round(x2[b, 0], 12)); k = K.index(round(x3[b, 0], 12))
        F[k*nk:(k+1)*nk, j*nj:(j+1)*nj] = np.array(a[b], dtype=float)
    return F


fv = sys.argv[1]; fh = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2].endswith('.bin') else None
fthin = float(sys.argv[-1]) if not sys.argv[-1].endswith('.bin') else 0.5
d = bc.read_binary(fv)
f1, f3, e = full(d, 'vet_f1'), full(d, 'vet_f3'), full(d, 'vet_e')
a1 = np.abs(f1).sum(axis=1)                     # (k, i)
p = f3.sum(axis=1)/np.where(a1 > 0, a1, 1.0)    # (k, i)
ff = a1.sum(axis=0)/np.maximum(e.sum(axis=(0, 1)), 1e-300)
thin = np.where(ff >= fthin*ff.max())[0]
odd = np.abs(p - p[::-1]).max(axis=0)/2.0
print('%s  t=%.6e  shells %d, thin %d (i %d-%d)' % (fv, d['time'], len(ff), len(thin),
      thin.min(), thin.max()))
print('  F3/|F1| odd: thin max %.3e median %.3e | thick(i<thin) max %.3e | thin mean <p> %.2e'
      % (odd[thin].max(), np.median(odd[thin]), odd[:thin.min()].max() if thin.min() > 0
         else 0.0, p[:, thin].mean()))
if fh:
    h = bc.read_binary(fh); rho = full(h, 'dens'); vz = full(h, 'velz')
    print('  v_phi rms: thin %.4e  mass-weighted all %.4e  max|v_phi| thin %.4e'
          % (np.sqrt((vz[:, :, thin]**2).mean()), np.sqrt((rho*vz**2).sum()/rho.sum()),
             np.abs(vz[:, :, thin]).max()))
