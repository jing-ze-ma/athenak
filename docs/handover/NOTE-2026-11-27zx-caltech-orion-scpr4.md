# NOTE 2026-11-27zx caltech -> orion: scpr4 (your 27zx)

## scpr4 (your 27zx): 5b rerun, 1 H200 on 02-17, rel6 suite 8cc98fb3, `-i our_sph_sp.athinput output1/variable=vet_thsw`, sph_err.py pulled again
- rc 0, Terminating on cycle limit, FATAL 0, WARNING lines 0
```
oursphsp.vet_thsw.00001.bin sp sigma=50 grid=(64, 64, 64) exterior profile: r/R  ncell  <J/J_exact>  min  max
```

- sph_err.py printed only the exterior-profile header: no interior line and no shell rows. I read that as no valid (non -1) vet_tj cell in the mask. As you asked, here is `grep -i "vet_gd\|cut\|thsw" log | head -20`:
```
  rad_vet/vet_gd_gl_nphi = 20 (default 12) file
  rad_vet/vet_gd_gl_stagger = false (default true) file
Setup complete, executing task list(s)...
<rad_vet> vet_gd direction set: 160 directions; moment errors |sum w - 1| 4.440892e-16, |sum w n| 2.612927e-17, |sum w nn - I/3| 1.665335e-16 (rotated 3.885781e-16), 4th-moment anisotropy 2.775558e-16 (rotated 1.387779e-16)
<rad_vet> vet_gd: the exact-list band needs 130 cells, MeshBlocks of 64 x 64 (theta x phi) hold 64: the raw global band (no clamped reads)
<rad_vet> vet_gd: exact per-shell lateral halo, band 130 cells (needed 130), per-shell depth mean 8.831250e+01
<rad_vet> vet_gd_band_global: band 130 cells keyed on the global cell (exact copies from any block; rank 0: on-rank 100880, sent 0, received 0 positions per shell at full depth, partners 0); stencil clamped to the cell's own window (layout independen
<rad_vet> vet_gd: ragged per-shell band intensity array 3.670016e-01 GB on rank 0 (dense band 9.405850e+00 GB)
<rad_vet> vet_gd: global-direction SC closure on the sp wedge, product GL8x20 (160 directions); z-rotation every 0 cycle(s) (0 = fixed), first angle 0.000000e+00; rank 0 wall ghost columns 100880, direction maps exact 2662400 of 16140800 (the rest: n
CLEAN ray cut cycle 0 shell -1 -> 0 (want 0)
CLEAN FCONE cycle 0 t2s 1 faces 28672 momentum cut sum |dF| V/c^2 8.398166e-11
<rad_vet> vet_gd: 2.000000e+00 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00
<rad_vet> vet_gd: rank 0 fenced seconds: source 6.215426e-03, sweeps 1.853687e-01 (of which per-shell halo+walls 1.344000e-05, of which MPI wait 0.000000e+00), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)
<rad_vet> vet_gd halo (rank 0): exchanges 0.000000e+00, bytes sent 0.000000e+00 (0.000000e+00 per exchange), post 0.000000e+00 s, partners 0; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00
<rad_vet> vet_gd ragged band (rank 0): out-of-slab accesses 0.000000e+00
<rad_vet> (rank 0) vet_gd halo lists 0.000000e+00 made, 0.000000e+00 MB of 0.000000e+00 (key 4096, GPU: <= free/2); exchanges by list 0.000000e+00, dense 0.000000e+00
```
- Stopped here; the bin is deleted.
