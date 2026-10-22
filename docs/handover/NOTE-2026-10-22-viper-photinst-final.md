# NOTE 2026-10-22 viper -> caltech: photosphere instability, viper FINAL (RSG E deck; work stopped, Caltech owns the dt-collapse test)

## 1. Seed-free static arms (he_seed 0, no MLT scaffold, radiation coupled): lateral rms of v_r at 0.99 / 0.995 / 1.0 R_ph [cm/s]
| grid | arith | kirchhoff | harmonic |
|---|---|---|---|
| cs N=32, 1-7.4 d | 1.3-2.4e3 / 1.3-3.5e2 / 42-190 | 1.3-39 / 1.3-23 / 0.4-2.3 | 0.24-4.0 / 0.03-0.55 / 0.04-0.12 |
| sp wedge 2.5 deg, 2-15 d | 1.7-3.5e3 / 2.8-5.1e2 / 63-470 | 0.18-2.7 / 0.36-1.5 / 0.04-1.35 | 0.08-3.5 / 0.12-0.56 / 0.07-0.26 |
- Either face option suppresses the bare-cliff symmetry breaking by 1e2-1e4 on BOTH grids; harmonic is lowest. A sub-cm/s residual
  remains that correlates weakly with cell geometry (corr ~0.1-0.17 with dOmega on cs, ~0.5-0.6 with sin(theta) on sp).
- The m=0 infall of the top layers (no scaffold) is the same in all arms (e.g. -2.0e5 cm/s at 0.99 R_ph, 7.4 d).
- The interior stays symmetric to round-off in every arm (and in pure-hydro cs_test, CPU).

## 2. RSG E N=32 restart arms from 40 d, at 57 d (harmonic's last dump 58.1 d)
| | base (arith) | kirchhoff | harmonic |
|---|---|---|---|
| corr(L 1.2-1.4 R_ph, dOmega) | +0.138 | +0.081 | +0.074 |
| corr(T_eff, dOmega) | +0.063 | +0.047 | +0.043 |
| corr(L, seed pattern) | -0.86 | -0.96 | -0.96 |
| cliff-shell e/rho sd (0.995 R_ph) | 0.148 | 0.013 | 0.012 |
| dt retries | 100 (from 56.1 d; dt 1000 -> 415-447 s) | 0 (dt 1000 s to 68.9 d) | 0 (dt 1000 s to 57.9 d) |
- The surface pattern is mostly the fixed he_seed function; its geometric part is small and lower with both face options.
- No dt collapse with harmonic at N=32 up to 58 d (your RSG E collapse was at 107 d, N=128, which these arms did not reach).
- Gates: kirchhoff and harmonic key-off bitwise vs f55be371 (50 cycles, hst); key-on C==C2 and R==C bitwise for both.

## 3. Earlier corrections stand
The interior "cube-locked" pattern was the fixed he_seed (57f45227). A symmetric star on cs stays symmetric to round-off; the cliff breaks
symmetry on sp too, so the cliff instability is not cs-specific (ee3efe63).

## 4. Not done
MI300A cs_test build segfaults at setup for iprob 13 (pure-hydro arms ran on CPU); zoom-grid key (local commit 620464da7) built but never run.
Files: /viper/ptmp2/jinma/photinst_1022 (scr/*.py; logs/static_face_cs.txt, static_face_sp.txt, imprint_*.txt; runs/*).
Also: accretor production run12 released on Raven 10-22 with fin51-r-fix ff1f52bc + implicit_face_opac = harmonic + time2_vet_col = extrap2.
