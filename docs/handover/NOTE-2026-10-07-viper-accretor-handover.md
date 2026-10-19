# NOTE 2026-10-07 viper -> DeltaAI: accretor hand-over and the env_rho_ph answer

Reply to TASK-2026-10-07-viper-accretor-handover.md.

## 1. What is in docs/handover/accretor-handover/

Verbatim copies of viper's `/viper/ptmp2/jinma/accretor_1006` small text files: DESIGN.md, IMPL.md, IMPL2.md,
CANDIDATES.md, PLASKETT.md, PENDING_CMDS.txt, msg*.txt, every input (bin/), the keys.txt of every arm
(keys/<tree>_<arm>.txt), the submit scripts (jobs/), the analysis/design scripts (scripts/ incl.
roche_stream.py and grid_design.py, tests/, page/, cand/) and the lit/ file list. **Start with INDEX.md**: run
history (every arm, job id, outcome, lesson), commits fe16e65a / 6588011b / ea04cd19 and binaries env6/env7,
science quantities, literature, limitations, next steps, and the paths of the big files left on viper.

Viper jobs 12116831 (q12_s1) and 12116832 (q12f_s1) were **cancelled while pending at the user's request
(10-07)**, after DeltaAI's copies had started. Viper will not run or change accretor work from now on.

## 2. How env_rho_ph = 0.3 was derived

**Recorded derivation (PLASKETT.md, "Density ratio (ESTIMATE, Mdot 1e-4)")**, the same formula as DESIGN.md
for RY Per:

    rho_stream = Mdot / (v A),   A = (c_s,don/Omega)^2 = (1.53 Rsun)^2,   v ~ 650 km/s
    rho_ph     = P_ph / c_ph^2,  P_ph = (2/3) g / kappa_e

Numbers (recomputed here, they reproduce the note): Mdot = 1e-4 Msun/yr = 6.30e21 g/s; A = 1.13e22 cm^2;
rho_stream = 8.6e-9 g/cc ("~9e-9"); g = G 16 Msun / (9.0 Rsun)^2 = 5.42e3 cgs, kappa_e = 0.34, P_ph = 1.06e4;
c_ph = 20 km/s -> rho_ph = 2.65e-9 g/cc; ratio 0.31 -> **env_rho_ph = 0.3**.

So, to the question:
- The cross-section is **(c_s/Omega)^2, not 2 pi (c_s/Omega)^2**, and treats c_s/Omega as the full width
  in both directions (no Gaussian factor).
- The speed is **v ~ 650 km/s, stated without derivation**. It is the RY Per number (DESIGN.md: "6e6 K for
  the full 657 km/s", the RY Per impact speed) carried over; the note calls the result the stream density
  "at impact". The Plaskett pgen gives other speeds (run.log of gpu2/q12_s45): |v| = 428 km/s rotating frame
  (500 inertial) at R_acc, and **135 km/s at r_out** (v_r -123.05, v_phi 56.08), where the code density unit
  rho_stream = 1 is actually imposed (peak of the phi Gaussian in the r_out window).
- **Physical input held fixed:** the notes record Mdot ~ 1e-4 as the QUOTED input (Wade+2026, PLASKETT.md
  table) and env_rho_ph as the derived ESTIMATE; IMPL2.md describes env_rho_ph as "the single ratio" that
  stands in for Mdot ("env_rho_ph = 20 corresponds to Mdot ~ 1e-6"). That is option (a). The notes do not
  record any decision about what to hold if the width changes.

RECONSTRUCTION (not in the notes; labelled as mine):
- Keeping the note's formula and only changing A to sigma_y sigma_z = 0.54 x 0.541 Rsun^2 raises the stream
  density 8.0x: env_rho_ph = 0.3/8.0 = **0.037** at Mdot 1e-4; holding 0.3 instead means Mdot = 1.25e-5.
- A consistent version (Gaussian peak, 2 pi sigma_y sigma_z, at the speed where the code imposes it, 135
  km/s) gives rho_ph/rho_peak = 0.40 for the old width (close to 0.3 only by coincidence of the 2 pi and
  the 650/135 speed factors) and **0.051** for 0.54 x 0.541. With 650 km/s and the Gaussian: 0.24.
- In the code the model is 2-D r-phi: the Gaussian is in phi only, the 4-cell theta band is uniform with
  reflecting faces, so at fixed code rho_stream the code Mdot drops only by the phi width ratio
  1.529/0.54 = 2.83, not by 8. Whatever sigma_z is assumed enters only through the physical-units mapping
  (env_rho_ph), not through the code.
- The default width uses the donor sound speed from t_don/mu_don (20.963 km/s), not env_cs_stream. The
  Ryu+2025 widths are at L1; the window is imposed at r_out = 0.85 d_L1 where the stream has accelerated to
  135 km/s, so its width there is not the L1 width (not computed on viper).

**Everything else in the setup that assumed the old width / density scale** (plaskett_env10..12f):
1. `stream_nsig = 3` window: half-width 3 sigma = 19.5 deg at r_out (old), 6.9 deg (new); outside it the
   outer BC is outflow-only. Resolution at r_out: r dphi = 0.0414 Rsun, sigma = 37 cells (old), 13 (new).
2. `env_amb_rho = 1e-6`, `env_amb_k = 3`, `env_sponge_rho = 10` (sponge on rho < 10 rho_amb), `rho_amb` =
   `hydro/dfloor` = 1e-7 ("1e-7 of the stream peak"): all in stream-peak units. The IMPL2.md ambient
   criteria (P_amb << P_ph, P_amb << stream ram, ambient mass << stream inflow) were checked against the
   old ratio; the Plaskett q10 failure interface had a 2e4 density contrast (dense stream layer vs ambient)
   that grows with the stream density.
3. `hydro/wb_rmax = 8.52` = R_acc - 1.5 d_pen with d_pen ~0.32 Rsun from ram-pressure balance at
   rho_s/rho_ph = 1/0.3 and 800 km/s (PLASKETT.md table: ram/surface pressure 5.3e3, penetration 0.32 Rsun).
   Under (a) the ram pressure rises ~8x and the penetration is deeper; under (b) unchanged.
4. Radial grid: fine zone 8.80..9.13 (dr 2.5-2.8e-3 = H_p/4) is set by H_p, not by the width; whether it
   covers the deeper penetration under (a) is unchecked.
5. dt: in every Plaskett run the limiting cell was hot ambient (c_ad ~390 km/s), dt 1.8e-6 -> 9.3e-7 as the
   stream arrives (q11_s45); not set by the stream width directly. `time/dt_min = 1e-7` collapse guard.
6. Stream inflow in code units (Min, dMin) and the ambient-mass-vs-inflow budget scale with the phi width;
   j_acc and torque/Mdot ratios do not.
7. Input comments still say "env_amb_rho ... (stream peak 1, rho_ph 20)" (RY Per) and carry the RY Per
   header.
