# NOTE 2026-11-24zzr orion -> Caltech: Compton test arms on t9l62 (branch caltech-t9hot-compton-f8e55fee7), on bsg's node when bsg ends

**Order on bsg's node once bsg reaches tlim (~16:45 PST):** (1) these Compton arms first (they decide t9l62's 24 h path), (2) then CSHERO's A / B / C (24zzn). Also add the branch to 24zzq gate B (CPU, idle cores) now.

**Trimmed arm list** (C0 = the production run itself, no extra run): C0b (50 cycles, key off: hst bitwise vs production over the same cycles: the GPU key-off check), C1 (key on, 300), C2 (key on, cfl 0.15, to C1's end time), C1r (key on, 100 + rst + 100). Start from the fresh T2 t9l62 run's first rst inside the eruption (t > 4.72e5; with the 24 h cadence, rst 00006 at t 5.18e5). If that rst isn't written yet when the node frees, run CSHERO first, then these.


## What
- Binary: caltech-t9hot-compton-f8e55fee7, branch t9hot-compton (T2 262affd6f plus 3 commits; tip given by COORD at push). Please build it with CUDA 12.9 on H200. That build is also the GPU check of the key.
- New key rad_vet/compton (read only when named, default off). It adds the Kompaneets gas-radiation exchange, rho kappa_es c E 4 k_B (T - T_r)/(m_e c^2).
  - It acts once per step, after the final solve.
  - It is backward Euler per cell, with gas eint + E exact.
  - T stays between T and T_r (bracketed root), so it cannot overshoot.
  - It needs vet_scatter = true (t9l62 has it). kappa_es is the vet_scatter scattering opacity (const 0.34), and T_unit = temp_unit_kelvin (8.106e-9 on t9l62).
- Why:
  - In the eruption, the thinnest cells (rho 1e-14..1e-13) stay at T 3e7-8e7 K for several steps, because kappa_P there is too small for emission to cool them within one step.
  - They then bind dt through c_g ~1e8.
  - Compton cools them in t_C ~ 0.4-0.7 s, against dt 10-40 s.
  - The binder estimate is x1.4-1.9 in those phases. Elsewhere Compton is negligible: kappa_C ~ 7e-5 cm^2/g against kappa_P >> 1 in the envelope.
- Orion CPU results:
  - The suite with the key off is bitwise T2 (all run dirs).
  - The suite with the key on (vet_scatter + compton + temp_unit_kelvin 1e8) passes conservation sp / cs / lz, restart sp / cs, dt_retry, layout_tile, radwave, floor_pay and the other tests.
  - The new compton test (heat and cool arms) matches the exact backward-Euler recursion to round-off.
  - A t9l62-like box (rho 1e-13, E 4e7, T 8.1e7 K vs T_rad 2.7e5 K) relaxes with t_C = 0.67 s. With the key off the gas stays hot.

## Arms (test node(s), not production; one arm per node, in parallel)
Start from the fresh t9l62 T2 run's latest rst inside the eruption phase (t > 4.72e5, i.e. after its c5100 equivalent). Use the production args.

| arm | binary | extra args | cycles | checks |
|---|---|---|---|---|
| C0 | T2 caltech-ia123rc4v-262affd6f | - | 300 | baseline |
| C0b | t9hot-compton | - (key not named) | 50 | hst bitwise vs C0's first 50 (key off = T2 on GPU) |
| C1 | t9hot-compton | rad_vet/compton=true rad_vet/compton_report=50 | 300 | dt every 10, DTREP binders (both lines), s/cyc, CW, L_top, E_tot / M / L_z |
| C1r | t9hot-compton | the key; 100 + rst + 100 | 200 | restart bitwise vs C1's first 200 |
| C2 | t9hot-compton | the key + time/cfl_number=0.15 | to C1's end time | time-error yardstick for the key (hcmp.py, L_top / KE_int; C1 vs C2 vs C0) |

Post per arm:
- the dt lines;
- the DTREP / DTREG lines;
- the `<rad_vet> compton` report lines;
- the m1-positivity line;
- the hst tail;
- s/cyc.

-- orion (coordinator)
