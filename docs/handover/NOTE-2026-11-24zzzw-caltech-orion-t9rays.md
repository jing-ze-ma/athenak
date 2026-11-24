# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
# NOTE 2026-11-24zzzw caltech -> orion: t9l62 ray-set test G1 / G2 (your 24zzzt)

- Preflight G2 (13x36, +1 cycle): 2026-11-24 13:35:39 rc=0 wall=57 s CW 0 FATAL 0. 
- Start t9l62.00015 (cycle 5279), T2 binary, C0 args plus your key list, 02-11 (7333073). Arms run: G1 and G2.

## G1
- 2026-11-24 13:42:15 rc=0 wall=395 s CW 0 FATAL 0; CW 0, FATAL 0, DT_RETRY 0; s/cyc 
- startup lines:
```
<rad_vet> vet_gd direction set: 126 directions; moment errors |sum w - 1| 6.66134e-16, |sum w n| 1.99493e-17, |sum w nn - I/3| 2.77556e-16 (rotated 1.11022e-16), 4th-moment anisotropy 2.498e-16 (rotated 1.11022e-16)
<rad_vet> vet_gd: exact per-shell lateral halo, band 57 cells (needed 57), per-shell depth mean 56.9688
<rad_vet> vet_gd: ragged per-shell band intensity array 17.4638 GB on rank 0 (dense band 61.1577 GB)
<rad_vet> vet_gd halo (rank 0): exchanges 7560, bytes sent 6.65277e+09 (879996 per exchange), post 0.126336 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0
```
- box stats (17x15x21 of gid 0; see assumptions at the end):
```
dump names ['cycle', 'dens', 'dt', 'eint', 'lx', 'nblk', 'time', 'velx', 'vely', 'velz', 'vet_d11', 'vet_d12', 'vet_d13', 'vet_d22', 'vet_d23', 'vet_d33', 'vet_e', 'vet_ex1', 'vet_f1', 'vet_f2', 'vet_f3', 'vet_fx1', 'vet
cycle t dt | cells rho<1e-11 | max|v| p99|v| (rho<1e-11) | n fast (|v|>5e7) | kF0/c median p90 [cm/s^2] | F0/(cE) median p90
5280 4.75224e+05 13.94 | 15 | 2.319e+08 2.174e+08 | 15 | 7.219e+05 2.651e+06 | 5.323e-02 1.944e-01
5281 4.75245e+05 20.99 | 15 | 1.683e+08 1.623e+08 | 15 | 7.002e+05 2.513e+06 | 5.184e-02 1.842e-01
5282 4.75267e+05 21.75 | 16 | 1.496e+08 1.473e+08 | 16 | 6.879e+05 2.253e+06 | 5.103e-02 1.652e-01
5283 4.75293e+05 25.9 | 16 | 1.536e+08 1.532e+08 | 14 | 7.205e+05 1.791e+06 | 5.282e-02 1.297e-01
5284 4.75318e+05 24.6 | 14 | 1.231e+08 1.227e+08 | 9 | 6.725e+05 1.909e+06 | 5.022e-02 1.383e-01
5285 4.75345e+05 27.12 | 14 | 1.198e+08 1.190e+08 | 9 | 8.524e+05 1.957e+06 | 6.125e-02 1.418e-01
5286 4.75373e+05 28.21 | 14 | 1.235e+08 1.234e+08 | 11 | 1.014e+06 2.041e+06 | 7.285e-02 1.494e-01
5287 4.75401e+05 28.57 | 13 | 1.431e+08 1.418e+08 | 10 | 1.325e+06 2.176e+06 | 9.583e-02 1.573e-01
ASSUMPTIONS: kappa 0.34; F0 = F - (4/3) E v (isotropic pressure; the exact form with D from vet_vet is in the npz for Orion). unit check F/E max = 1.468e+09 (F energy flux if > 1)
```
- user.hst (t dt L_top/L_in):
```
4.752244e+05 13.94 2.55543
4.752454e+05 20.99 2.55536
4.752671e+05 21.75 2.55505
4.752930e+05 25.9 2.55477
4.753176e+05 24.6 2.55453
4.753447e+05 27.12 2.55428
4.753729e+05 28.21 2.55401
4.754015e+05 28.57 2.55372
```

## G2
- 2026-11-24 13:43:33 rc=0 wall=62 s CW 0 FATAL 0; CW 0, FATAL 0, DT_RETRY 0; s/cyc 
- startup lines:
```
<rad_vet> vet_gd direction set: 468 directions; moment errors |sum w - 1| 1.77636e-15, |sum w n| 2.22804e-17, |sum w nn - I/3| 2.22045e-15 (rotated 2.05391e-15), 4th-moment anisotropy 3.85803e-15 (rotated 1.47105e-15)
<rad_vet> vet_gd: exact per-shell lateral halo, band 57 cells (needed 57), per-shell depth mean 56.9688
```
- box stats (17x15x21 of gid 0; see assumptions at the end):
```
dump names None
cycle t dt | cells rho<1e-11 | max|v| p99|v| (rho<1e-11) | n fast (|v|>5e7) | kF0/c median p90 [cm/s^2] | F0/(cE) median p90
Traceback (most recent call last):
  File "/resnick/groups/carnegie_poc/jingze/t9rays_1124/ext.py", line 52, in <module>
    print('ASSUMPTIONS: kappa 0.34; F0 = F - (4/3) E v (isotropic pressure; the exact form with D from vet_vet is in the npz for Orion). unit check F/E max = %.3e (F energy flux if > 1)'%fn)
NameError: name 'fn' is not defined
srun: error: hpc-sm-02-11: task 0: Exited with exit code 1
```
- user.hst (t dt L_top/L_in):
```
```

- Boxes (.npz per cycle, with 2 cells each side): caltech_results/T9RAYS/{G1,G2}/. The local indices are read as 0-based, gid 0 = first block in the dump (lx stored). The full dumps are deleted after extraction.
