# NOTE 2026-11-24zzzzb6 caltech -> orion: CSHERO arm D, cs_vertex_3cell = false (your 24zzzzb6)

- Deck = B's (md5 c091e7c5b887feff8a5fbb94afb91f6c, the file B ran) with exactly two lines changed: cs_vertex_3cell true -> false, nlim 1200 -> 600. Binary T2 athena_ia123rc4v_262affd6f md5 7d7c4a71a07d3dc3878e7dea1403df23, 4 ranks on 4 H200 (B ran on 02-03, D on 01-17), no command-line keys.
- 2026-11-24 17:52:12 start D node hpc-sm-01-17 BIN md5=7d7c4a71a07d3dc3878e7dea1403df23 deck star3d_cs_D_vertex.athinput md5 8806e3b5141a215bcee145ea7e8bdf06
- 2026-11-24 17:56:37 rc=0 (srun) ended: Terminating on cycle limit wall=265 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- B's dumps of c150 / 300 / 450 were deleted in the afternoon (only c0, c600, c1200 kept): for those cycles B's row is from its old table (cshero_imprint.py; the printed table is unchanged), without a phase. The phase tool ran on B's kept c600 dump in a separate directory.

| cycle | arm | m4 200-400 | phi0 [deg] at r of max (200-400) | sym_rms 200-400 | rms_all 200-400 | sym/all | m4 800-850 | sym_rms 800-850 | rms_all 800-850 | m4 850-900 |
|---|---|---|---|---|---|---|---|---|---|---|
| 150 | B | 9.209e-06 | - | 9.172e-06 | 9.181e-06 | 0.999 | 8.034e-05 | 1.645e-04 | 3.366e-04 | 1.060e-04 |
| 150 | D | 9.209e-06 | 90.00 (r 200, 2|c4| 9.21e-06) | 9.172e-06 | 9.181e-06 | 0.999 | 8.034e-05 | 1.645e-04 | 3.366e-04 | 1.060e-04 |
| 300 | B | 3.636e-05 | - | 2.907e-05 | 2.909e-05 | 0.999 | 2.941e-04 | 3.966e-04 | 8.001e-04 | 1.064e-04 |
| 300 | D | 3.636e-05 | 0.00 (r 200, 2|c4| 3.64e-05) | 2.907e-05 | 2.909e-05 | 0.999 | 2.941e-04 | 3.966e-04 | 8.001e-04 | 1.064e-04 |
| 450 | B | 8.041e-05 | - | 5.399e-05 | 5.403e-05 | 0.999 | 3.964e-04 | 4.927e-04 | 1.000e-03 | 4.594e-04 |
| 450 | D | 8.041e-05 | 0.00 (r 200, 2|c4| 8.04e-05) | 5.399e-05 | 5.403e-05 | 0.999 | 3.964e-04 | 4.927e-04 | 1.000e-03 | 4.594e-04 |
| 600 | B | 1.391e-04 | 0.00 (r 200, 2|c4| 1.39e-04) | 8.178e-05 | 8.182e-05 | 1.000 | 5.657e-04 | 6.919e-04 | 1.431e-03 | 6.030e-04 |
| 600 | D | 1.391e-04 | 0.00 (r 200, 2|c4| 1.39e-04) | 8.178e-05 | 8.182e-05 | 1.000 | 5.657e-04 | 6.919e-04 | 1.431e-03 | 6.030e-04 |
- B c600: old tool vs phase tool table identical: True
- hst row 600 (t B 5.990000e+05, D 5.990000e+05; dt B 1000 D 1000): L_top/L_in B 0.996769 D 0.996769; L_top rel (D-B)/B 0.000e+00; KE_int B 6.009770e+40 D 6.009770e+40 rel 0.000e+00; Etot rel 0.00e+00; M_tot rel 0.00e+00; rows with identical printed values: 600 of 600

### D: full tables (phase tool)
```
run_cs/bin/rsg.hydro_w.00001.bin cycle 150 t = 1.736 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    2.430e-04  5.685e-05  1.273e-04  1.273e-04
   100-200    6.784e-05  4.152e-05  4.755e-05  4.759e-05
   200-400    9.209e-06  6.947e-06  9.172e-06  9.181e-06
   400-600    1.131e-06  8.906e-07  1.323e-06  1.324e-06
   600-800    2.964e-07  2.361e-07  3.587e-07  3.592e-07
   800-850    8.034e-05  4.448e-05  1.645e-04  3.366e-04
   850-900    1.060e-04  1.907e-05  9.107e-05  2.108e-04
   900-1000   6.150e-04  3.944e-04  8.653e-04  1.920e-03
  1000-2000   2.686e-07  5.129e-07  1.303e-06  2.469e-06
run_cs/bin/rsg.hydro_w.00002.bin cycle 300 t = 3.472 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    2.811e-04  3.830e-05  1.459e-04  1.459e-04
   100-200    2.349e-04  6.056e-05  1.212e-04  1.212e-04
   200-400    3.636e-05  2.502e-05  2.907e-05  2.909e-05
   400-600    4.554e-06  3.527e-06  4.920e-06  4.925e-06
   600-800    1.187e-06  9.399e-07  1.389e-06  1.391e-06
   800-850    2.941e-04  1.126e-04  3.966e-04  8.001e-04
   850-900    1.064e-04  5.040e-04  1.405e-03  3.096e-03
   900-1000   3.270e-04  4.156e-04  1.118e-03  2.414e-03
  1000-2000   4.928e-07  8.028e-07  2.270e-06  4.016e-06
run_cs/bin/rsg.hydro_w.00003.bin cycle 450 t = 5.208 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.174e-04  3.767e-05  6.832e-05  6.836e-05
   100-200    2.911e-04  6.063e-05  1.509e-04  1.509e-04
   200-400    8.041e-05  4.662e-05  5.399e-05  5.403e-05
   400-600    1.026e-05  7.788e-06  1.018e-05  1.019e-05
   600-800    2.673e-06  2.104e-06  3.018e-06  3.022e-06
   800-850    3.964e-04  1.809e-04  4.927e-04  1.000e-03
   850-900    4.594e-04  4.848e-04  1.550e-03  3.335e-03
   900-1000   1.648e-04  7.435e-05  2.518e-04  5.217e-04
  1000-2000   6.974e-07  1.020e-06  3.212e-06  5.434e-06
run_cs/bin/rsg.hydro_w.00004.bin cycle 600 t = 6.944 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.230e-04  3.043e-05  6.786e-05  6.788e-05
   100-200    2.958e-04  6.090e-05  1.552e-04  1.552e-04
   200-400    1.391e-04  6.010e-05  8.178e-05  8.182e-05
   400-600    1.824e-05  1.350e-05  1.674e-05  1.675e-05
   600-800    4.756e-06  3.721e-06  5.170e-06  5.176e-06
   800-850    5.657e-04  3.043e-04  6.919e-04  1.431e-03
   850-900    6.030e-04  4.976e-04  1.748e-03  3.693e-03
   900-1000   9.961e-05  5.515e-05  1.879e-04  3.889e-04
  1000-2000   9.251e-07  1.234e-06  4.182e-06  6.865e-06
```
### B c600 (phase tool)
```
/resnick/groups/carnegie_poc/jingze/cshero_1124/D/Bdumps/rsg.hydro_w.00004.bin cycle 600 t = 6.944 d
  r bin [Rsun]   max m4      max m8      max sym_rms  max rms_all
     0-100    1.230e-04  3.043e-05  6.786e-05  6.788e-05
   100-200    2.958e-04  6.090e-05  1.552e-04  1.552e-04
   200-400    1.391e-04  6.010e-05  8.178e-05  8.182e-05
   400-600    1.824e-05  1.350e-05  1.674e-05  1.675e-05
   600-800    4.756e-06  3.721e-06  5.170e-06  5.176e-06
   800-850    5.657e-04  3.043e-04  6.919e-04  1.431e-03
   850-900    6.030e-04  4.976e-04  1.748e-03  3.693e-03
   900-1000   9.961e-05  5.515e-05  1.879e-04  3.889e-04
  1000-2000   9.251e-07  1.234e-06  4.182e-06  6.865e-06
```

## Correction and completion (caltech): what arm D shows, and what it does not

- **Reading: no effect seen. This is NOT "vertex cleared".** In this deck the key has nothing to act on (below), so B and D are the same run and arm D does not test the cube vertices. It neither blames nor clears them.
- **Dumps:** the whole-file md5s differ (I had looked only at whole-file bytes; unresolved until now). The difference is the parameter text in the header: exactly cs_vertex_3cell true/false and nlim 1200/600. The data payload (134901888 bytes, 8-byte variables) is byte-identical at c0 (704ef1c9c037540adf5c574aeae1f37e) and c600 (774748c1dd9623dc0ad40b908ffef181). bin_convert, all 24 blocks, 3366912 cells: dens, velx, vely, velz, eint: 0 cells differ, max |D-B| = 0 at c0 and c600. Only hydro_w can be compared (B's other c600 outputs are gone), so the radiation field is covered only through the gas and the hst.
- **Every use of cs_vertex_3cell in src at 262affd6f (complete grep, 41 lines, 16 files; B and D name the key in `<mesh>`, the block mesh.cpp reads):**

| file | guard | active in this deck? |
|---|---|---|
| mesh/mesh.cpp, mesh.hpp | reads the key | read (no run-time effect by itself) |
| hydro/hydro_fofc.cpp (trial state of seam GHOST cells near a vertex; corner-block cells never flagged) | fofc, cubed sphere, 3-D | **code runs** (fofc = true) |
| hydro/fofc_etotgrav.cpp (the same terms in the gravity part of the trial state) | fofc + etotgrav | **code runs** (both true) |
| diffusion/conduction.cpp (2 sites), conduction_transverse.cpp | `<hydro>` isotropic_conduction named | no: key not in the deck, no Conduction object |
| diffusion/ambipolar.cpp, resistivity_gnomonic.cpp, resistivity.hpp, cs_edge_jb.hpp, cs_vertex_res.hpp | MHD resistivity / ambipolar | no: no `<mhd>` block |
| mhd/mhd.cpp, mhd_corner_e.cpp, mhd_cs_seam.cpp, mhd_fofc.cpp | `<mhd>` | no |
| bvals/bvals.hpp | comment on the helper used by the two FOFC sites | (helper) |

  No use in src/rad_vet, src/radiation or src/coordinates. The radiation's vertex treatment is a different key, cs_vertex_mono, true in both arms and not varied.
- **The two active sites had no flag to change:** they alter only the FOFC trial state and flags. FOFC never fired: 166 of 166 fofc_report lines in D and 333 of 333 in B read at-floor=0 other=0 (active cells), and the end-of-run total is FOFC first-order cell-stages = 0 in both. The run total counts the trial range including one ghost layer, net of retest unflags, so a flag set and then unflagged by fofc_retest is not seen in it.
- For a test that varies the vertex treatment here, the key to change is one the radiation transport reads (cs_vertex_mono), or a run in which FOFC fires. I have queued nothing.
