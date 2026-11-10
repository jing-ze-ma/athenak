# TASK 2026-11-10 viper -> ORION: raycheap18 -- the dt-consistent half-range blend (implicit_blend_dtc) on the sp suite + cs Part 2

Answers NOTE raycheap17 (e01dcddea) and raycheap17b-rootcause (8a1d985bc). CPU-small, as raycheap17. Times PST.

## 0. What viper found and changed (branch fork/raycheap18-1110 = d981e435c, one commit on dea3e3732)
Same mechanism as your 17b section 1: the stored x1/lateral face flux is the blend (1-w) th (F^n - ...) + w F_hr with F^n
the previous BLENDED flux, so a fixed w relaxes the face memory by the fraction w PER STEP (rate w/dt). Consequence
measured on viper (pulse_u SP T, profiles): the fine-dt solutions carry a period ~3 cell oscillation in i 42-59 (E 0.17
<-> 0.25) that GROWS as dt shrinks (BE damping of the central-flux grid modes falls with dt): not an offset, a dt-damped
grid mode. Steady state of the old blend: weight on F_hr = w (1 + a)/(w + a), a = c dt chi_f -> depends on dt.

Fix (opt-in, no new memory array): `implicit_blend_dtc = kappa > 0` (read only when named; absent = bitwise) replaces the
per-step weight of every half-range face (M1HrFace and M1ThinFace; x1 kernel and the lateral / cs face body) by
  wt = w a/(a + 1 - w),   a = c_hat dt (chi_f + kappa/dx) = (c dt/dx)(tau_f + kappa).
With kappa -> 0 the steady face flux is EXACTLY the convex blend (1 - w) F_diff + w F_hr for every dt, and the
relaxation rate c chi_f w/(1-w) is dt-free (a fixed ODE as dt -> 0); wt -> w when c dt chi >> 1 (production thick /
large-step regime nearly unchanged), wt(0) = 0, wt(1) = 1, wt <= w (M-matrix kept). Cost: one division per face.
This is the algebraic equivalent of your design (a) at steady state, without carrying F_c through t2f1 / t2i / kk / restarts.
Use kappa = 1e-6 (kappa = 1 makes the vacuum rate c/dx dominate: pre-asymptotic, see table).
Design (b) (half-range at the stage time): the extrap2 machinery already has it: `time2_vet_x2_fields = 13` (9 + bit 4 =
vgd_hr extrapolated to t^{n+1} like the tensor; production 9 leaves vgd_hr at t^n).

## 1. viper measured (CPU login, binary /viper/ptmp2/jinma/builds/bin/athena_none_cpu_d981e435_rc18, work dir
/viper/ptmp2/jinma/rc18_1110/w, scripts = raycheap17 st17.py + your pulse_u vimp drop + EXTRA/TAG hooks)
| arm (pulse_u, E) | d per pair | p |
|---|---|---|
| SP T (reproduces yours) | 3.01e-2 3.58e-2 3.79e-2 3.97e-2 3.24e-2 | -0.25 -0.08 -0.07 0.29 |
| SP + dtc 1 T | 1.49e-2 1.16e-2 7.30e-3 4.30e-3 2.00e-3 | 0.36 0.67 0.76 1.10 |
| SP + dtc 1e-6 T (l to 64) | 8.83e-3 5.12e-3 2.32e-3 9.15e-4 2.77e-4 1.19e-4 | 0.79 1.14 1.34 1.72 1.22 |
| **SP + dtc 1e-6 + x2_fields 13 T** | 9.26e-3 5.09e-3 2.11e-3 7.24e-4 1.74e-4 | **0.86 1.27 1.55 2.06** (F1 ... 1.99) |
| SP + dtc 1e-6 C | 6.33e-2 6.12e-2 5.73e-2 4.67e-2 | 0.05 0.09 0.30 -- SPACE still fails (see 3) |
Gates: keys absent bitwise vs dea3e3732 (SP T_1/2/4 tab files, 17 digits, identical); S09+gd with dtc on bitwise
(T_1, T_2: the closure beam does not take the half-range path).
(x2_fields 13 alone, without dtc, does not cure the flat T: dtc 1 + x2 13 gave 0.51 0.72 0.68 1.02.)

## 2. Build (Orion)
`git fetch fork raycheap18-1110` (d981e435c), the raycheap14 CPU build line with that commit. Keys-absent gate: SP pulse_u
T_1 and fs_u C_32 bitwise vs athena_dea3e3732.

## 3. Runs (st17.py, your drops; new arms via EXTRA as in your hooks)
FIX = `implicit_blend_dtc=1.0e-6 time2_vet_x2_fields=13`; FIXN = FIX + `vet_col_ncore=64`.
(a) pulse_u: SP-FIX T l 1..64, C, X (levels to 512); SP-FIXN C, X; references SP, S09+gd (yours). Pass: T p >= 1.8 at the
    finest pairs; C/X as 09-25.
(b) SP + vet_col_ncore 64 (C, X) on pulse_u and fs_u (already in your 17 note for C; add X).
(c) Part 1b: pulse_s, pulse_edd_u, lat_th_atm, lat_ph_atm, mvp_u with SP, SP-FIXN, S09 (C, X, T). Report every variable.
(d) SPACE bisection on pulse_u X with FIXN (the C failure above is not the dt mechanism): one-off arms
    implicit_thin_corr=none; implicit_hr_face=cell; implicit_thin_switch=none + implicit_blend_flo=0 implicit_blend_fhi=1e-6;
    vet_gd_rotate_every=0 vet_gd_rotate_dt=0 vet_rotate_avg=1 vet_rotate_avg_ema=0; implicit_thin_js_face=none (or its
    off value); implicit_thin_hcface=mean. Name the first that restores p ~2 in X.
(e) RELEASE the held cs Part 2 of raycheap17 with CSP and CSP-FIXN arms (same EXTRA; the key also acts on the cs
    lateral faces through rad_m1_latface_body.hpp).
Answer as NOTE-2026-11-1x-orion-raycheap18*.md on bsg-files-1009 with the eval lines; fast subset (a)+(d) first.
