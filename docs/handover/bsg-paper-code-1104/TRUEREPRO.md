# BSG TRUE reproduction of Ma, Bildsten & Jiang 2026 (arXiv 2609.28656): setup prepared, production NOT submitted (10-08)

Goal (user 10-08): a reproduction in which the ONLY difference from the paper is the radiation transport
(our implicit grey M1 + VET closure instead of Athena++ implicit discrete ordinates, Jiang 2021).

- Work dir: `/viper/ptmp2/jinma/bsg_truerepro_1008/`. All paths below are relative to it.
- Code: branch `bsg-truerepro-1008`, local commit `54bf6697` on `03590601` (= fork/rt-integration).
  - Worktree: `wt/`.
  - Not pushed, not merged. Only `src/pgen/he_star_m1.cpp` changed (+234/-16).

## 0. Summary

- **Matchable differences.** The finished "reproduction run" (repro_sf_4n -> viper repro -> repro_sf_cont, 57.4 d)
  differs from the paper in these settings, which we can match:
  - the well-balanced (WB) scheme and its reference force;
  - the inner hydro BC;
  - the outer hydro BC;
  - the flux surgery at both walls;
  - the 1e-2 seed;
  - the M1 initial flux above r_photo;
  - the radial grid;
  - the scattering split in the closure source.

  All of these are now set to the paper. The energy form of gravity (etotgrav) also differs but is kept: it is
  minor and the pgen needs it. Every other setting is the same, unknown (questions for Linhao), or intrinsic to
  the RT (section 1).
- **The paper's Fe peak and photosphere are its evolved 3-D state, not its IC** (section 2).
  - Paper IC (1-D): Fe peak ~40-41 Rsun, photosphere ~51 Rsun.
  - Paper at 40 d (3-D): Fe peak ~35 Rsun, photosphere 43 Rsun. The envelope contracted by 6-8 Rsun and released
    ~5e44 erg.
  - Ours: IC 41.35 and 53.0 Rsun; at 30-57 d 40.1 and 51.9 Rsun. It hardly contracted (dE_grav -4.7e43 erg,
    dU +4.3e43 erg).
- **1-D columns: without WB the column does NOT hold HSE as well as with WB** (section 4).
  - The no-WB arms drift secularly. At 11.6 d (3.6 Fe-peak thermal times):
    - rho is -7.8 % at 20.5 Rsun and -12 % at 35 Rsun;
    - T at the Fe peak is -7.6 %;
    - the Fe-peak isotherm moves inward 0.77 Rsun;
    - the wall-cell rms v_r is 35x the WB arm's.
  - WB arm: 41.351 -> 41.310 Rsun.
  - The paper's scheme (no WB) therefore carries a drift that WB suppressed in our runs. It goes in the direction
    of the paper's contraction.
  - In 1-D the zero-contact BC and our old wall give the same result to <= 0.5 % in rho. The WB switch dominates.
- **The zero-contact inner BC holds the mass to round-off.** The integrated mass through the inner face is 3.2e-15
  of the column mass over 11.6 d.
- **Thin-atmosphere signature (coordinator request): PRESENT in the finished BSG reproduction** (section 6).
  - Above 1.05 R_ph the transported face flux is saturated (|F| >= 0.99 cE in 71-100 % of faces).
  - 38-69 % of those faces point inward.
  - Lateral face fluxes are 50-400x the shell-mean radial flux.
  - This is the AG Car central-flux pattern. It is a candidate for the photosphere offset and **must be fixed
    before the true reproduction is launched.**
- **GPU smoke: PASSED.** apudev job **12132658**, 10 cycles clean, same binary and input.
- **Cost.** 19-28 h on 8 viper apu nodes, or 24-37 h on 4 Raven nodes (section 7).
- **Submit commands: not executed** (section 8).

## 1. Difference table: paper vs the finished reproduction run

Sources:
- The paper text: `bsg_1001/2609.28656.txt`, sections 2.2.1-2.2.5.
- The restart-embedded input of the last viper rst: `bsg_viper_1006/repro/rst/bsg3d.00051.rst`, extracted to
  `diff/rst00051_input.txt`.
- XKEYS of every link: `rad_m1/time_scheme=be time/cfl_number=0.3 rad_m1/force_reference_work=split` + path
  overrides (`bsg_viper_1006/repro/run.cfg`, `merged_repro/raven_cont_run.cfg`).

| # | setting | paper | finished run | mark | true-repro input |
|---|---|---|---|---|---|
| 1 | RT method | implicit discrete ordinates, 120 angles, backward Euler, modified-HLLE face flux (alpha 5), tol 1e-6 | implicit grey M1, vet_gd T5wihc closure, be, tol 1e-8 | DIFFERENT (RT-intrinsic, keep) | unchanged |
| 2 | domain | 20-80 Rsun, theta pi/3-2pi/3, phi 0-pi/3, periodic angles | same | SAME | same |
| 3 | radial grid | base 192 (text "196"; 2560 blocks of 16^3 needs 192) + SMR x2 on 35-55 -> dr 0.15625 / 0.3125 Rsun | nx1 267 smooth plateau, 0.153 on 34-55.9, 0.306 outside | DIFFERENT (can match) | **nx1 256: 0.15625 on 35.0-55.16, 0.3125 outside** (`grid/`) |
| 4 | angular grid | 256^2 in 35-55, 128^2 elsewhere | 256^2 everywhere | DIFFERENT (code: rad_m1 on sp has no SMR); ours finer, keep | unchanged |
| 5 | IC hydro | radiative-only MESA profile (Linhao file) + isothermal atmosphere to the floor, v = 0 | same file, our reconstructed r_photo 53.03 Rsun | SAME as far as known; r_photo UNKNOWN (Q1) | same |
| 6 | IC radiation | two-zone intensities: E = aT^4, F = L/4pi r^2 below r_photo, **F = 0 above** | F = L/4pi r^2 everywhere | DIFFERENT (can match) | **`ic/ic_ma2026_f0.txt`** |
| 7 | seed | none | 1e-2 T seed, k 2-4, 16 modes, 37.5-45.2 Rsun (added 10-06) | DIFFERENT (can match) | **he_seed 0** |
| 8 | EOS | ideal mu 0.62, gamma 5/3 (+ radiation) | same | SAME | same |
| 9 | gravity | point mass 20 Msun, Athena++ source term, no WB | point mass in the **x1 WB polytropic pair** + reference force a_ref = kappa F_IC/c (wb_arad, split work) + WB guard/clamp | DIFFERENT (can match) | **he_gravity = pointmass, wellbalance_dynamic = false, force_reference = none** |
| 10 | energy form of gravity | Athena++: energy source from mass fluxes | etotgrav (rho Phi in E, conservative) | DIFFERENT (numerics, minor; pgen needs etotgrav) | unchanged |
| 11 | opacities | TOPS X 0.7 Z 0.008 (Linhao tables), edges clamped; kappa_a = kR - ke, ke 0.34, kP | Linhao tables padded (= clamp); flux kR, exchange kP | SAME (moment level) | same |
| 12 | scattering in the closure source | Jiang 2021 eq. 6, uncapped | capped eps = min(kP/kR, 1) | DIFFERENT (RT, can match; effect <= 6e-6, VET_SCATTER.md) | **vet_scatter = true** |
| 13 | Riemann solver | HLLC | HLLC | SAME | same |
| 14 | reconstruction | not stated (Athena++ PLM?) | PLM (Mignone non-uniform) on WB deviations | UNKNOWN (Q4); WB part removed | PLM, no WB |
| 15 | hydro integrator | not stated (vl2?) | rk2 | UNKNOWN (Q4) | rk2 |
| 16 | CFL | 0.3 hydro (0.025 "radiation") | 0.3 with rad_signal_speed (dt up to ~1.6x smaller) | UNKNOWN (Q10); keep for stability | same |
| 17 | FOFC | not stated | on (pgen default) | UNKNOWN (Q11); keep (guard, fires mostly in floor gas) | same |
| 18 | floors | "density floor", value not given | dfloor 1e-20, pfloor 1e-12, tfloor 5000 K, vceil 1e8 | UNKNOWN (Q3) | same |
| 19 | inner BC hydro | ghosts = UNSCALED IC rho, P; ghost 1: v_r solved for HLLC S_M = 0 after reconstruction, v_t = 0; other ghosts v = 0 | closed mirror wall (wall_closed_ix1), ghosts = IC scaled by the edge cell, v mirrored, plus face-flux surgery (he_wall_zero_flux) | DIFFERENT (can match) | **he_bc_inner = zerocontact, he_wall_closed = false, he_wall_zero_flux = false** |
| 20 | inner BC radiation | ghost two-zone intensities from T_IC(r_g), F = L/4pi r_g^2; DO modified-HLLE face flux -> thick limit ~ (F_g + F_a)/2 + (c/4 tau_c)(E_g - E_a), bath term <~ 4 % | M1 face flux fixed = L/4pi r_in^2 | DIFFERENT (RT-intrinsic; same L in steady state; paper's has a weak T-bath restoring term) | unchanged |
| 21 | outer BC hydro | copy rho, P, v_t; v_r copied if > 0 else 0 | noinflow: ghosts = scaled IC, inflow mirrored, plus inflow-flux removal | DIFFERENT (can match) | **he_bc_outer = copy, no surgery** |
| 22 | outer BC radiation | vacuum: outgoing copied, incoming 0 | Marshak, E_bath 0, VET surface q | DIFFERENT (RT-intrinsic; moment form of vacuum) | unchanged |
| 23 | MLT scaffold / IC re-balance | none / none | mlt_flux_frozen false / he_ic_balance false | SAME | same |
| 24 | sponge, reduced c | none | none, chat/c 1 | SAME | same |
| 25 | M1 guards | n/a | implicit_resid_fatal 1e-2, pos_gas/pos_floor, realisable_coupling (fix A), tsolve_opac, gas_newton_switch 10, thin_freeze | RT-intrinsic, keep | same |
| 26 | nghost | 2 (assumed for xorder 2) | 3 | UNKNOWN, minor (zero-contact uses ghost 1; v = 0 in 2-3) | same |
| 27 | run length | steady after 30 d, to 57.4 d | 57.4 d | SAME | tlim 4.96e6 s |

Notes:
- Row 25: resid_fatal only stops a diverged solve. The positivity guards act only where E or the gas energy would
  go negative (floor and thin cells). Neither changes the physics where it does not fire.
- Row 16: changes dt by at most a few tens of %.

## 2. Fe peak and photosphere: IC or evolved state?

`fepeak/fepeak.py`: kappa max for r > 25 Rsun; tau = 1 integrated from the top.

| profile | Fe peak | r(tau = 1) |
|---|---|---|
| Linhao file, own (MESA/OPAL) kappa | 41.10 Rsun (178 kK, kappa 1.43) | 48.32 (the file ends at 52.49 with no photosphere) |
| Linhao file, Ma TOPS kR | 41.35 Rsun (174 kK, kappa 1.59) | 48.32 |
| our repro IC (`ic_ma2026.txt`) | 41.35 | 52.89 (r_photo built at 53.03) |
| paper 1-D (text, Fig. 3 right) | ~40 | ~51 |
| paper 3-D, t = 40 d (Figs. 3, 5) | ~35 | 43 (42-44) |
| our finished run, 30-57.4 d (`cmp57.json`) | 40.12 | 51.95 (51.3-52.7) |

- 35/43 Rsun is the paper's EVOLVED state: inward by ~6 Rsun at the Fe peak and ~8 Rsun at the photosphere.
- That is consistent with dE_grav ~ -5e44 erg. A rough estimate: 0.01 Msun moved from 45 to 38 Rsun gives ~3e45
  erg.
- Ours barely moved: Fe peak 41.35 -> 40.1 Rsun.
  - dE_grav(25 d) = -4.7e43 erg.
  - dU = +4.3e43 erg: the envelope GAINED internal energy.
  - L_top/L_in = 0.955.

## 3. Implementation (opt-in, default off)

`wt/src/pgen/he_star_m1.cpp`, commit 54bf6697.

- **`problem/he_gravity = wb (default) | pointmass`.**
  - Requires etotgrav = true, wellbalance_dynamic = false and force_reference = none. Refuses he_gm_column and
    he_ic_balance.
  - Momentum source, the Athena++ PointMass form:
    `d(rho v_r) = -dt rho GM coord_src1/x1v`, with `coord_src1 = 1.5 (r+^2 - r-^2)/(r+^3 - r-^3)`.
  - Energy: etotgrav carries rho Phi, so no energy source is needed.
  - No WB effective potential and no reference force are set up.
- **`problem/he_bc_inner = zerocontact`.**
  - All inner ghosts hold the UNSCALED IC (rho, eint), v = 0.
  - Exception, ghost is-1: v_t = 0 and v_r = vg, the root of S_M(vg) = 0. S_M is the HLLC contact speed at the
    wall face, from:
    - the left state reconstructed from cells (is-2, is-1, is) and the right state from (is-1, is, is+1);
    - both made with the same `Hydro::PLM_nonuniform` as the x1 sweep;
    - then hllc_hyd.hpp steps 2-5 (`HsZcAm`).
  - Root find: bracket expansion, then Illinois regula falsi to |S_M| <= 1e-14 c_s. Fallback with no bracket:
    mirror v.
  - S_M = 0 makes the HLLC mass, transverse-momentum and energy fluxes vanish; the momentum flux is p*.
  - Requires spherical polar, an ideal gas, hllc, reconstruct_x1 plm, no WB, he_ic_balance false, nghost >= 2.
    It defaults he_wall_closed to false and refuses true.
- **`problem/he_bc_outer = copy`.** Every ghost copies the edge rho, eint and v_t. v_r is copied if > 0, else 0.
  M1 ghost fill as for the other outer BCs.
- **Gates:**
  - Key off is **bitwise** against base 03590601 on CPU (`col/bit_base` vs `col/bit_new`, 30 cycles of col_wb):
    hst identical, 10 dumps identical.
  - cpplint: the same 3 messages as the base file.
  - Lines <= 90 columns.
  - No new compiler warnings.
- **Binaries:**
  - CPU: `/viper/ptmp2/jinma/builds/bin/athena_he_cpu_54bf6697_truerepro`.
  - GPU: `/viper/ptmp2/jinma/builds/bin/athena_he_gpu72_54bf6697_truerepro`, md5
    978765ae88b0634b3bb876430a0869ef (copies in `smoke/bin/`, `prod/bin/`).

## 4. 1-D column tests

Setup (`col/`):
- Login node, 1 rank each, nice.
- Repro N267 grid, 4 x 4 angular cells.
- RT keys of the column check + implicit_realisable_coupling.
- t = 0 -> 1e6 s = 11.6 d = 3.6 thermal times at the Fe peak (U(r > 41.35)/L = 3.2 d).

Arms:
- **wb**: the finished run's setting.
- **zc**: true-repro.
- **nowb_wall**: no WB with our old wall/noinflow.

Analysis: `col/col_ana.py`; output in `col/col_ana_now.txt`. Shell-mean change since t = 0, rho / T:

| at 11.57 d | wb | zc | nowb_wall |
|---|---|---|---|
| 20.5 Rsun (wall region) | +0.9 % / +0.26 % | -7.8 % / +5.3 % | -7.3 % / +5.0 % |
| 25 Rsun | -0.68 % / -0.04 % | -0.81 % / -1.6 % | -1.1 % / -1.6 % |
| 35 Rsun | -2.0 % / -0.06 % | -12.3 % / -2.4 % | -12.5 % / -2.4 % |
| 41.35 (Fe peak) | -3.0 % / -0.35 % | -4.2 % / -7.6 % | -4.3 % / -7.6 % |
| 45 Rsun | +0.7 % / -0.7 % | -34 % / -16 % | -34 % / -16 % |
| r(T = 1.739e5 K), IC 41.351 | 41.310 | 40.578 | 40.576 |
| rms v_r of the wall cells (last) | 14 m/s | 490 m/s | 498 m/s |
| L_top/L_in (last 10 %) | 0.942 | 0.895 | 0.896 |
| Etot - Etot0 vs int (L_in - L_top) dt | 1.007e39 / 1.013e39 | 1.525e39 / 1.541e39 | 1.522e39 / 1.528e39 |
| int Mdot_bot dt / M | 0 (closed) | 3.2e-15 | 0 (closed) |
| dominant Mr_int period (t > 1e5 s) | 1.74 d | 1.74 d | 1.74 d |

What this shows:
- **HSE without WB: holds worse.** The drift is secular: near the wall rho -0.7 %/d with T rising, and the
  wall-cell rms v_r is 35x larger.
  - The Fe-peak isotherm moves inward 0.77 Rsun (WB: 0.04).
  - The initial settling is the same in both arms (~1.2-1.6 % at 25 Rsun at 0.23 d): the IC is ~1 % out of HSE
    with mu 0.62 + point mass.
  - So the difference is the scheme. Athena++ without WB at the paper's resolution should carry a similar drift.
    Its direction (inward, energy released) matches the paper's contraction.
- **Inner BC.** Zero-contact holds the mass to round-off. In 1-D it matches our old wall to <= 0.5 % in rho. Its
  effect on the 3-D transient (bouncing radial flows) is a 3-D question.
- **Radial mode.**
  - Same dominant period, 1.74 d, in all arms. Resolution is coarse: a 10.4 d window, so neighbouring bins are
    2.1 / 1.45 d.
  - The no-WB arms have 1.2x the rms radial momentum.
  - The 3-D 3.91 d mode does not appear in 4 x 4 columns.

## 5. 3-D input and GPU smoke

- Input: `run3d/bsg3d_truerepro.athinput`.
  - Built from the finished run's master input (`repro_final/bsg3d_repro.athinput`).
  - Rows 3, 6, 7, 9, 12, 19 and 21 of the table set to the paper; seed keys removed; viper paths.
  - Extra output7 `m1_face` every 0.5 d for the atmosphere check.
- IC: `ic/ic_ma2026_f0.txt`, md5 246a666e57abe5052828c271d4ff3ca4.
- CPU check of the N256 grid with zc: `col/smk_zc256`, 20 cycles, rc 0, dx1 at the wall 0.3125 Rsun.
- **Smoke `smoke/smoke.sh`, apudev job 12132658** (vipa1327, 1 node x 2 MI300A, 8 blocks/GPU):
  - Env exports HSA_XNACK=1 and HSA_NO_SCRATCH_RECLAIM=1. The guarded apu duplicate 12132659 was cancelled while
    pending.
  - **(A) true repro, same binary and input, t = 0, 10 cycles.**
    - rc 0, 0 FATAL/NaN.
    - Picard mean 4.6, max 7, NON-CONVERGED 0. Non-fatal Newton-fallback notes at i 168-250.
    - dt 87.3 s, the same as the finished run's start.
    - 16-25 s/cycle on 2 GPUs.
    - vet_scatter echoed on.
    - user.hst at t = 874 s: Mdot_bot 4.9e10 g/s (round-off), wall rms v_r 89 m/s, L_top/L_in 1.00.
  - **(B)** the finished run's binary (d6164a6a) from viper rst 00051 (t = 48.2 d), 3 cycles, writing m1_face for
    section 6. rc 0, 40-75 s/cycle on 2 GPUs.
  - Monitor (`gpumem.12132658.log`): GPU memory-allocated column 91 % on both devices. This is the monitor's
    VRAM-allocated fraction (%), not GB. Production uses 8 nodes.
- Smoke-only overrides: nlim, ndiag, output7/dt. The 10-cycle smoke of the chain script itself (section 8) is the
  formal gate before the production chain.

## 6. Thin-atmosphere check on the finished BSG reproduction (coordinator request)

**(a) Cell-centred dumped flux** (`atmchk/atm_sig.py` -> `atm_sig.txt`):
- Data: merged_repro dumps 61-117, every 8th, 30-57 d; R_ph = 51.945 Rsun.
- m1_f1 is lab-frame, face mean + v(E+P), capped at cE.

| r/R_ph | 0.95 | 1.0 | 1.05 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 |
|---|---|---|---|---|---|---|---|---|
| shell L/L | 0.954 | 0.853 | 0.153 | 0.014 | -0.128 | -0.267 | -0.207 | -0.127 |
| inward cells | 0 | 0.043 | 0.41 | 0.49 | 0.56 | 0.68 | 0.56 | 0.60 |
| \|F\| >= 0.99 cE | 0 | 0.13 | 0.93 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| lateral E rms | 0.29 | 0.74 | 0.076 | 0.065 | 0.058 | 0.062 | 0.056 | 0.074 |

**(b) Transported (comoving) face flux** (`atmchk/face_sig.py`, `face_sig_o.py` -> `face_sig_B*.txt`):
- Data: smoke B's m1_face at t = 4.1613e6 s, with E from its m1 dump.

| r/R_ph | 0.9 | 0.95 | 1.0 | 1.05 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 |
|---|---|---|---|---|---|---|---|---|---|
| 4 pi r^2 <F_x1>/L | 0.957 | 0.955 | 0.954 | 0.960 | 1.02 | 6.05 | 13.6 | 10.6 | 2.21 |
| inward faces | 0 | 0 | 0.027 | 0.38 | 0.50 | 0.57 | 0.69 | 0.56 | 0.61 |
| \|F_x1\| >= 0.99 cE | 0 | 0 | 0.053 | 0.71 | 0.98 | 1.00 | 1.00 | 1.00 | 1.00 |
| rms(F_x2, F_x3)/\|<F_x1>\| | 0.05 | 0.12 | 10 | 48 | 232 | 190 | 52 | 101 | 395 |

- The outer top face (fx1o) carries 0.946 L. hst L_top: 0.955.

**Verdict: present.**
- From 1.05 R_ph to the top (1.54 R_ph) the transported field is saturated beams, about half or more pointing
  inward.
- The lateral face fluxes are ~10^2 x the radial mean.
- The comoving shell L is 6-14 L at 1.2-1.4 R_ph. The top face still gives the right L. I did not check whether
  the lab advective term v(E+P) cancels the excess.
- Below 1.0 R_ph the field is clean.
- The radiation force on the atmosphere gas comes from these fluxes, so the state above the photosphere is wrong.
- This is a candidate for the photosphere offset (51.9 vs 43 Rsun) and the weak top variability.
- **It must be fixed before the true reproduction.** Candidate: the sp blend flux that does not show it in the AG
  Car tests. It is not on this branch. Merge it, rebuild, re-smoke this input, and re-run `face_sig.py` on the
  smoke's m1_face.

## 7. Cost estimate

Basis: the finished run, `repro_final/cost.json`.
- 57200 cycles to 57.4 d, mean dt 86.7 s.
- Raven 4 nodes x 4 A100: 1.59-1.60 s/cycle.
- viper 8 apu nodes x 2 MI300A: 1.22 s/cycle.

| platform | finished run | true repro (grid x 0.96; dt x 1.0-1.5 if paper-like flows, KE 8x, shorten dt) |
|---|---|---|
| viper 8 nodes | 19.4 h, 155 node-h | **18.6-28 h, 150-225 node-h** (2 x 24 h links) |
| Raven 4 nodes | 25.4 h, 102 node-h | **24.4-37 h, 98-147 node-h** |

- The WB walk is removed, so the per-cycle cost may drop a little (not measured). The 2-GPU smokes are
  memory-bound and are not a scaling basis.
- Raven needs a CUDA build of the branch. It is local only, so a push or an rsync of the tree is needed (user
  decision). Then build in an interactive job.

## 8. Submit commands (NOT run)

Before anything: fix section 6, rebuild, and update `prod/run.cfg` BIN/BIN_MD5. Then smoke the chain script with
the same binary, input and command line plus NCYC:

```
mkdir -p /viper/ptmp2/jinma/bsg_truerepro_1008/prod_smoke && cp /viper/ptmp2/jinma/bsg_truerepro_1008/prod/run.cfg /viper/ptmp2/jinma/bsg_truerepro_1008/prod_smoke/
sbatch -J bsg_trs -p apudev --nodes=1 --time=00:14:00 \
  -o /viper/ptmp2/jinma/bsg_truerepro_1008/prod_smoke/link.%j.out -e /viper/ptmp2/jinma/bsg_truerepro_1008/prod_smoke/link.%j.err \
  /viper/ptmp2/jinma/bsg_truerepro_1008/prod/truerepro.sh /viper/ptmp2/jinma/bsg_truerepro_1008/prod_smoke 10
```

Production on viper: 8 nodes, 2 x 24 h chained afterany.

```
D=/viper/ptmp2/jinma/bsg_truerepro_1008/prod
J1=$(sbatch --parsable -J bsg_tr --nodes=8 --time=24:00:00 -o $D/link.%j.out -e $D/link.%j.err $D/truerepro.sh $D) && \
J2=$(sbatch --parsable -J bsg_tr --nodes=8 --time=24:00:00 -d afterany:$J1 -o $D/link.%j.out -e $D/link.%j.err $D/truerepro.sh $D) && echo $J1 $J2
```

`prod/truerepro.sh`:
- A fresh start from IN when `rst/` is empty; otherwise a restart from the newest rst with output last_time reset
  (rst_info.py).
- Guard, CANCEL/DONE/STOP, md5 check, NaN/FATAL -> STOP.
- Exports HSA_XNACK=1 and HSA_NO_SCRATCH_RECLAIM=1.
- Syntax-checked (`bash -n`), not run.

## 9. Open questions for Linhao

Earlier list (MA2026_REPRO_IC.md sect. 5), still open:
- Q1: r_photo, T_photo and rho_photo of the IC, and the footnote-5 adjustment.
- Q2: atmosphere formula: P total or gas, local g.
- Q3: the floors.
- Q4: reconstruction and hydro integrator.
- Q5: inner ghosts: unscaled? zero-contact only in layer 1? which ghost T?
- Q6: F = 0 above r_photo.
- Q7: sponge, velocity ceiling, reduced c.
- Q8: HSE with mu 0.62.
- Q9: the SMR region.

New:
- Q10: does the hydro dt include radiation pressure in the signal speed? What does "CFL 0.025 for radiation" mean
  in the implicit scheme?
- Q11: FOFC or another positivity fallback? NGHOST?
- Q12: "196 base cells" vs 2560 x 16^3 blocks: is the base 192, i.e. SMR exactly 35-55 Rsun at dr 0.15625?
- Q13: the zero-contact root: tolerance; once per step or per RK stage; with that stage's reconstructed states?
- Q14: was the 3-D contraction (Fe peak 41 -> 35, photosphere 51 -> 43 Rsun, dE_grav -5e44 erg) checked against
  resolution, or with a hydrostatic test without convection? Our no-WB columns drift the same way.
- Q15: the radiation inner ghost: F at the ghost centre? alpha = 5 in the modified HLLE?
- Q16: the light curve: taken at the outer boundary or at tau = 1? (For the 16x variability gap.)

## 10. Files

- `diff/rst00051_input.txt`.
- `fepeak/fepeak.py`.
- `ic/ic_ma2026_f0.txt`.
- `grid/faces_N256_ra35.0_w0.45.txt`.
- `col/`: col_{wb,zc,nowb_wall}.athinput, run.sh, col_ana.py, col_ana_now.txt; arms wb/ zc/ nowb_wall/; smk_*;
  bit_base/ and bit_new/.
- `run3d/bsg3d_truerepro.athinput`.
- `smoke/`: smoke.sh, A/, B/, smoke.12132658.out, gpumem.12132658.log.
- `atmchk/`: atm_sig.py/.txt, face_sig.py, face_sig_o.py, face_sig_B*.txt.
- `prod/`: truerepro.sh, run.cfg, bin/.

## 11. Blend switch (10-08/09, merge2-1008 agent): prepared, NOT submitted

- **Input.** `run3d/bsg3d_truerepro.athinput` now sets `implicit_flux = blend` and `implicit_blend = tau_f`.
  The central copy is `run3d/bsg3d_truerepro_pre_blend.athinput`.
- **Binary.** `prod/run.cfg` BIN = `prod/bin/athena_he_gpu72_ee03881b_merge2`, md5 1c7365f2f4b7681f93d07c06c5361cf9.
  - This is merge2-1008 = sp-blend 4622f424 + truerepro 54bf6697 on rt-integration cc7e4948. It was never pushed.
  - **SUPERSEDED** by rc-1009/accel-1009. Before production, rebuild at accel-1009 + truerepro, update BIN/BIN_MD5,
    and re-smoke.
  - Dry-run merge of 54bf6697 onto accel-1009 (c72ada3d): one conflict in `src/pgen/he_star_m1.cpp`, in the
    gravity source. It is additive: keep the fixbundle F6 tally block, then the `if (hs_pm_) {...} else {` wrapper
    around hs_grav.
- **Chain smoke.** apudev 12134859, `prod/truerepro.sh prod_smoke2 10`, same input and binary.
  - Fresh t = 0: rc 0, 0 FATAL.
  - Picard mean 5.4, max 7, NON-CONVERGED 0.
  - dt 87.37 s; ~11 s/cycle on 2 GPUs.
- **Thin-atmosphere check: NOT DONE.**
  - A blend restart of the finished central repro rst 00051 (`smoke_blend/face_blend.sh`, apudev 12134860, without
    `implicit_recon_lag=step`) died on its first solve with FATAL `Picard solve DIVERGED` (resid 1.39e5). This
    matches AG Car (central -> blend restarts diverge).
  - The 10-cycle fresh smoke does not reach an output7 m1_face dump (dt 0.5 d).
  - Check production instead: run `python3 atmchk/face_sig.py prod/bin/bsg3d.m1_face.00001.bin <m1 dump at the
    same t>` on the first dumps.
- **Raven variant** (4 nodes x 4 A100, 16 ranks, one block per GPU), staged in /raven/ptmp/jinma/bsg_truerepro_1008:
  - `prod/truerepro_raven.sh`: port of truerepro.sh, with `-p gpu`, gres a100:4, KOKKOS_MAP_DEVICE_ID_BY=mpi_rank
    and `files/rst_info.py`.
  - `prod/run.cfg` and `prod_smoke/run.cfg`.
  - Input `files/bsg3d_truerepro_raven.athinput` (md5 f0e8d602), plus the IC and both tables with md5 equal to viper.
  - Binary `bin/athena_he_a100_ee03881b_merge2` (md5 dbc01d630161b1f1ff67487ed1594301; full CUDA build, interactive
    job 31011551). It is superseded the same way as the viper binary.
  - The Raven smoke was NOT submitted (the permission check refused it).
- **Submit commands** (only after the accel-1009 + truerepro rebuild and a re-smoke):

viper smoke, then production (8 nodes, 2 x 24 h):
```
D=/viper/ptmp2/jinma/bsg_truerepro_1008
mkdir -p $D/prod_smoke3 && cp $D/prod/run.cfg $D/prod_smoke3/
sbatch -J bsg_trs -p apudev --nodes=1 --time=00:14:00 -o $D/prod_smoke3/link.%j.out -e $D/prod_smoke3/link.%j.err $D/prod/truerepro.sh $D/prod_smoke3 10
P=$D/prod
J1=$(sbatch --parsable -J bsg_tr --nodes=8 --time=24:00:00 -o $P/link.%j.out -e $P/link.%j.err $P/truerepro.sh $P) && \
J2=$(sbatch --parsable -J bsg_tr --nodes=8 --time=24:00:00 -d afterany:$J1 -o $P/link.%j.out -e $P/link.%j.err $P/truerepro.sh $P) && echo $J1 $J2
```
Raven (ssh -S ~/.ssh/raven.sock raven01i.mpcdf.mpg.de) smoke, then production (4 nodes, 2 x 24 h):
```
D=/raven/ptmp/jinma/bsg_truerepro_1008
sbatch -J bsg_trs -p gpu --nodes=4 --time=00:15:00 -o $D/prod_smoke/link.%j.out -e $D/prod_smoke/link.%j.err $D/prod/truerepro_raven.sh $D/prod_smoke 10
P=$D/prod
J1=$(sbatch --parsable -J bsg_tr --nodes=4 --time=24:00:00 -o $P/link.%j.out -e $P/link.%j.err $P/truerepro_raven.sh $P) && \
J2=$(sbatch --parsable -J bsg_tr --nodes=4 --time=24:00:00 -d afterany:$J1 -o $P/link.%j.out -e $P/link.%j.err $P/truerepro_raven.sh $P) && echo $J1 $J2
```
- **Open: merged-binary wall time on the old repro restart.**
  - The keys-off merge2 binary ran the finished-repro rst 00051 (8 blocks/GPU, 90 % device memory) at 120-330
    s/cycle. 03590601 ran it at 35-70 s/cycle (jobs 12134854, 12135858).
  - The run.log arithmetic is identical. During the slow cycles the GPU sits idle for minutes.
  - The bisect job 12137385 was cancelled before it started (re-scope). Recheck with accel-1009. It may be memory
    oversubscription at 1 node; the true repro on 1 node ran at ~11 s/cycle.
