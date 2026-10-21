# TASK 2026-10-21 caltech -> raven: test vet_gd for the 3-D accretor (USER)

USER asks: test the vet_gd closure for the accretor.

- The current v12lzc_be_gm template (5de0635e) runs `closure = vet_col` with NO vet_gd keys. Its comment says
  "vet_gd cannot run on the 4-cell theta band (design sect. 4)" — that was the 2-D band. The production mesh is now 3-D:
  nx2 = 58 (+-20.39 deg theta wedge), nx3 = 512, nx1 = 420.
- Please check: (1) whether vet_gd supports this theta wedge with walls at both theta faces (sp mesh, not cs); (2) if not,
  what it needs; (3) if yes, a short accretor test with vet_gd (RSG E / HG keys as the reference set: twin, rotate, symref,
  seg, bc_diff, twin_fuse/twin_det) vs the vet_col baseline: cost per cycle, dt, redos/non-converged, energy budget, and how
  much the radiation field (F_r, E, f-factors) differs above the photosphere and in the stream.
- Then a template with vet_gd (md5) if it's a GO.
- Caltech status: accretor FRESH on 077dda23 + 5de0635e (vet_col) runs as 5223354 (stream on since 07:05 PDT, ~0.6 s/cycle,
  NOTE 5e85c2ab). It keeps running until you send a vet_gd template; then USER decides whether to restart fresh.
