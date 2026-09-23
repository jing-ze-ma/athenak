# Resolution and cost plan for the cubed-sphere deep hot Jupiter (prod4 and hyd4)

2026-09-22. This note is analysis only; no code was changed. Every number cites a file, a job
or a formula. Abbreviations: "prod3 dump" = `bench/cs_mhd_prod3/bin/dhj.mhd_w_bcc.00142.bin`
(rot 283.28). Measurements on that dump used `docs/handover/scripts/dhjcs.py` with the EOS
table `bench/cs_ens/analysis/eos_table.txt`. The scripts `res.py`, `nh.py` and `jet.py` are in
the session scratchpad (not kept). "Timing jobs" = 11943312 and 11943313 (see section 3).

## 0. The grid being assessed

From `inputs/production/deep_hot_jupiter_cs_prod4.athinput`. `deep_hot_jupiter_cs_hyd4.athinput`
has the same `<mesh>` and `<meshblock>` blocks.

| item | value |
|---|---|
| radius | r = 9.44e9 .. 2.0556e10 cm (r/a_p 1 .. 2.18) |
| radial cells | nx1 = 128, polynomial stretch c1..c4 = -0.068392, -2.191487, 2.464818, -1.366698 |
| horizontal cells | 6 panels x 32 x 32 (gnomonic equiangular); dxi = 90/32 = 2.81 deg; r dxi = 4.63e8 cm at a_p |
| MeshBlocks | 128 x 16 x 16, so 24 blocks (2 x 2 per panel), 786,432 cells |
| time | rk2, cfl 0.3, tlim 8.64e7 s = 283 rotations (P = 3.05e5 s) |
| MHD | eos resistivity, max_eta 5e12, no STS |
| RT | ck, ck_pcut_bar 10, spherical tm sweep, pseudo-spherical beam |

Radial cells, from the map `r = r0 + (r1-r0) u(xi)`: dr_min = 4.57e7 cm at r/a_p 1.171 and
dr_max = 1.86e8 cm at the top (`res.py`).

**Machine correction.** A viper node has **2 APUs**, not 4: `sinfo -p apu` reports `gpu:2`, 96 CPUs,
220 GB per node (325 nodes), Sockets=2; apudev likewise, MaxNodes=1. "1 node" = 2 GPUs below.

## 1. Physics length scales vs the grid

All values in this section are from the developed state in the prod3 dump.

### 1.1 Radial: pressure scale height vs dr

H_p = -(d ln p/dr)^-1 is computed per column; cells/H = H_p/dr (`res.py`, `nh.py`).

| p [bar] | median r/a_p | median H_p [km] | median dr [km] | cells/H_p median | p05 | p01 |
|---|---|---|---|---|---|---|
| 1e-6 (observable limb) | 1.732 | 9468 | 1195 | 7.8 | 2.7 | 1.6 |
| 1e-4 | 1.404 | 5066 | 746 | 6.8 | 2.6 | 2.1 |
| 1e-3 | 1.298 | 3794 | 586 | 6.5 | 2.2 | 2.1 |
| 2e-3 (H2 front) | 1.275 | 2730 | 552 | 5.0 | 2.2 | 2.1 |
| 2e-2 (H_p minimum) | 1.234 | 1452 | 499 | **2.9** | 2.5 | 2.4 |
| 0.1 | 1.209 | 1611 | 473 | 3.4 | 3.2 | 3.2 |
| 1 | 1.164 | 2040 | 458 | 4.5 | 4.4 | 4.3 |
| 10 (ck cut) | 1.108 | 2575 | 512 | 5.0 | 5.0 | 4.9 |

Over all cells with p > 1e-6 bar: median 5.0, p05 2.8, p01 2.3 cells per H_p (for H_rho:
6.0 / 2.7 / 2.1).

- **The stated requirement (10 cells/H below 1e-6 bar, memory `dhj-grid-resolution-design`) is
  not met.** The prod4 header's "median 10.1 / p05 7.3" belongs to the nx1 = 234 design.
- **The map refines the wrong radius:** fitted to the initial H_rho, it refines r/a_p 1.171
  (memory `radial-stretch-refit`); the developed H_p minimum is at r/a_p 1.23 (dr 495 km).
- **Cells needed:** the median column spans 19.3 H_p from the base (234 bar) to 1e-6 bar, plus
  4.4 above (`nh.py`). 10 cells per median H_p needs **191 cells below 1e-6 bar** (296 for the
  p05 H_p); the grid has 100 there and 28 above, where H_p ~ 11,600 km and dr ~ 1,480 km.

**Ck cut and photosphere.** At rot 283 the Rosseland cut is icut = 24-25 in every column
(median r_cut = 1.089e10 cm). The 10 bar level would give icut = 16 (`tests_ck_sph/REAL_RAYS.md`).
So the whole ck region is at 4.5-5 cells/H_p or better, except the 2e-3..0.1 bar band at 2.9-3.4.
The notes do not record the pressure of the thermal photosphere (tau = 2/3). I did not measure it.

**Convergence evidence (radial).** 40-cycle hydrostatic test: nx1 = 128 has 5x the spurious KE of
nx1 = 234 (5.94e32 vs 5.63e31, second order; memory `radial-stretch-refit`). 0.25-rot A/B uniform
200 vs poly 120 agreed within variability (memory `radial-grid-stretch`). The H2 front thickness
(2.07e8 cm) moved 1.02x for a 1.30x dr change, so it is physical (`dhj-grid-resolution-design`).
**No developed-state (tens of rotations) radial convergence test exists**; section 5 proposes one.

### 1.2 Horizontal: jet, deformation radius, vertices

**Cell size.** r dxi = 2.81 deg = 4.63e8 cm at a_p, and 6.0e8 cm at the 1e-3 bar radius. The
gnomonic line element is h = r dxi (1+X^2) sqrt(1+Y^2)/(1+X^2+Y^2), with X = tan xi and
Y = tan eta. The narrowest cell width is 0.724 r dxi, at the panel-edge midpoints and vertices.
The largest |cos| of the angle between grid lines is 0.475 (`res.py`).

**Deformation radius.** N^2 = g(d ln p/dr / Gamma1 - d ln rho/dr), with Gamma1 from the EOS
table, g = 942 (a_p/r)^2 and Omega = 2.06e-5. Medians per isobar (`res.py`):

| p [bar] | N [s^-1] | N H [cm/s] | L_eq = sqrt(NH r/2 Omega) | cells | L_D(30 deg) = NH/f | cells |
|---|---|---|---|---|---|---|
| 1e-4 | 5.6e-4 | 2.85e5 | 41 deg | 14.7 | 60 deg | 21 |
| 1e-3 | 6.8e-4 | 2.58e5 | 41 deg | 14.6 | 59 deg | 21 |
| 1e-2 | 8.6e-4 | 1.28e5 | 30 deg | 10.5 | 30 deg | 11 |
| 0.1, 1 | median N^2 <= 0 (neutral/convective) | - | - | - | - | - |
| 10 | 2.5e-5 | 6.4e3 | 7 deg | 2.5 | 1.7 deg | 0.6 |

The deformation radius is resolved by 10-21 cells in the stratified upper atmosphere. At
10 bar it is 0.6-2.5 cells, i.e. unresolved, but N there is 30x smaller and the flow is weak
(|u| < 0.2 km/s below 0.1 bar, memory `dhj-jet-shallow-westward-deep`).

**Equatorial jet width.** Time-mean zonal-mean u over 4 dumps at rot 168-174. sp is
`sp_mhd_nopole` (sinh theta stretch 3, so d(theta) = 0.84 deg at the equator, 3.3x finer than
cs meridionally). cs is `cs_mhd_prod3`. Source: `jet.py`.

| p [bar] | sp u_eq | sp FWHM | cs u_eq | cs FWHM | cs cells per FWHM |
|---|---|---|---|---|---|
| 1e-4 | 1.12 km/s | 10.1 deg | 0.92 km/s | 6.3 deg | 2.2 |
| 1e-3 | 0.96 | 11.3 | 1.33 | 11.3 | 4.0 |
| 1e-2 | 0.79 | 10.7 | 1.23 | 19.5 | 6.9 |

sp resolves the core with ~12 cells and finds it ~10 deg wide at every level; cs matches at
1e-3 bar, is narrower at 1e-4 bar (2.2 cells: grid scale) and twice as wide at 1e-2 bar. Not yet
a convergence result: sp also lacks rad_implicit_x1, rad_cap_ang, rt_use_cons (memory
`cs-vs-sp-dhj-rot50`). (The old half-strength cs-hllc upper jet was the solver: ausm+up on the same
grid matched sp, `dhj-jet-shallow-westward-deep`.) **Reading:** 32^2 gives 2-4 cells across the
core at p <= 1e-3 bar, marginal; 48^2 gives 3.3-6, 64^2 4.5-8.

**Terminator.** The pseudo-spherical beam deposits in a twilight ring down to mu0 = -0.45,
i.e. 27 deg beyond the terminator (`tests_ck_sph/REAL_RAYS.md`). Its dayside partner,
0 < mu0 < 0.45, is another 27 deg, so the ring spans ~9.5 cells on each side. The dominant
error in the ring is horizontal inhomogeneity along grazing rays: the beam deposits 0.5-0.99 of
the true power. That is a model limit, which refinement does not remove.

**Shocks and day-night flow.** Horizontal Mach v_h/c_s (c_s = sqrt(Gamma1 p/rho)), median / p99:
1e-4 bar 1.38 / 4.98 (max 5.4); 1e-3 bar 1.46 / 3.67; 1e-2 bar 0.78 / 1.93; 0.1 bar 0.16 median.

So the upper day-night flow is supersonic and its shocks are captured over 2-3 cells, not
resolved. Refinement sharpens them and moves them; no study in the notes measures what that
does to heat transport.

### 1.3 MHD: magnetic diffusion scale

Where eta is at the cap the resistive length is eta/U = 5e12/2e5 = 2.5e7 cm (U = 2 km/s at 1e-3
bar): 0.04 of a horizontal cell (Rm_cell = U h/eta = 2e5 x 6e8/5e12 = 24) and ~0.5 dr (radial
Rm_cell = |v_r| dr/eta = 0.1-4 for |v_r| 1e4-4e5 cm/s in the dump).

Where eta is below the cap (the ionised deep layers) the resistive length is smaller still.
The field is therefore limited by numerical resistivity horizontally at any affordable grid.
Evidence: the dump ME is 17-26 % below the face ME, i.e. the field is at grid scale (memory
`cs-vs-sp-dhj-rot50`). Horizontal refinement lowers the numerical resistivity but will not
converge Rm. MHD convergence has to be judged on integrated quantities: the zonal-KE brake
(31 %, memory `sp-hydro-vs-mhd-comparison`) and ME.

## 2. Timestep: which limit binds

dt is `cfl * min(dt1, dt2, dt3)`, a straight minimum over directions (memory `dt-binding-direction`).

- **Measured.** prod3 `dhj.mhd.hst` harmonic-mean dt by rotation band: 6.0 s (rot 0-1),
  10.4 s (1-5), 16.6 s (5-20), 20.4 s (20-50), 20.1 s (50-283). Minimum over rot 20-283:
  12.4 s.
- **Reconstructed.** A CFL-0.3 reconstruction on the prod3 dump gives **dt1 = 20.4 s
  (radial)**, which matches the measured 20.1 s (`res.py`).
  - The binding cell is at i = 45, r/a_p 1.257, 7.8e-3 bar: a downdraft with v_r = -4.1e5 cm/s
    and c_s = 3.6e5 cm/s at the H_p minimum.
  - The horizontal limit at the 0.724 vertex width is **dt_h = 62.6 s**, i.e. 3.1x slack. The
    worst horizontal cell is at 8e-6 bar, v_h 1.5e6 cm/s.
  - vA is at most 4.0e5 cm/s and only 3.9e4 cm/s at the binding cell, so the field does not
    set dt. Hydro and MHD share the same CFL dt.
- **Horizontal refinement is free in dt up to ~3x.** 48^2 gives dt_h ~ 42 s and 64^2 gives
  ~31 s, both above 20.1 s. Confirmed at t = 0 by the timing jobs: dt = 28.977 s for 32^2,
  48^2 and 64^2 alike (`run_*_{A,B,C,D}/run.log`). Caveat: a finer grid may develop faster
  small-scale downdrafts, which would lower dt1. Not measured.
- **Radial refinement scales dt as dr at the binding cell.** nx1 = 192 gives dt(t=0) =
  19.318 s, against 28.977 s at nx1 = 128: exactly 1.500x (`run_*_E/run.log`). The steady
  dt is therefore about 13.4 s.
- **MHD Ohmic limit.** dt_ohm = cfl (1/6) dr^2 / eta_cap (memory `dt-binding-direction`,
  exact to 7 digits).
  - Per the prod4 header, the cap is reached near r/R_p 1.24 (dr ~ 5.1e7 cm). With
    max_eta 1e13 that gave 12.9 s, so with 5e12 it gives 25.8 s at nx1 = 128.
  - At nx1 = 192 it gives **11.5 s < 13.4 s**, so the Ohmic term binds MHD. The fix is to
    scale max_eta with dr^2 (2.2e12) or turn on use_rkg_sts.
  - At the finest cell (dr 4.57e7) with 5e12 it would be 20.9 s, level with the CFL today.

## 3. Throughput: measured on the GPU

**Developed-state anchor, prod4 form at rot 283** (vipa1327, 300 cycles, 2 GPUs,
`tests_rank_imbalance/README.md`): 12/12 split **12.87 cycles/s** (mean of 6) = 1.012e7
zone-cycles/s; 10/14 split 14.73; tm arm of `bench/swform_0922` 14.35 (`tests_ck_sweep_form`).

**Scaling with GPU count**, prod3-restart config (`tests_mpi_gate_0922/README.md` B5). Per-GPU
zone-cycles/s:

| GPUs | nodes | blocks per GPU | zone-cycles/s per GPU | efficiency |
|---|---|---|---|---|
| 1 | 1 | 24 | 4.97e6 | 1.00 |
| 2 | 1 | 12 | 3.96e6 | 0.80 |
| 4 | 2 | 6 | 2.53e6 | 0.51 |

Going off-node costs 1.56x per GPU at this size. At 1.92M cells (08-26, memory
`meshblock-decomposition-gpu`) the numbers were 91 % at 2 GPUs and 69 % at 4 GPUs / 2 nodes.

**New timing jobs (this note):** **11943312** (arms A B C D E) and **11943313** (E D C B A), both
apudev node vipa1327; prod4 binary (md5 f95130b2, `bench/dhj_res_0922/athena`), prod4 input from
scratch, 150 cycles per arm, 2 ranks / 2 GPUs, outputs off. Code's own `zone-cycles/cpu_second`
(`bench/dhj_res_0922/run_<job>_<arm>/run.log`); memory = peak `rocm-smi` VRAM used (`mem.log`).

| arm | grid (nx1 x per-panel^2) | blocks x size | zones | zone-cycles/s, run 1 / run 2 | mean cycles/s | mem/GPU |
|---|---|---|---|---|---|---|
| A (prod4) | 128 x 32^2 | 24 x 128.16.16 | 786,432 | 1.024e7 / 1.067e7 | 13.29 | 6.8-8.9 GB |
| D | 128 x 48^2 | 24 x 128.24.24 | 1,769,472 | 1.411e7 / 1.402e7 | 7.95 | 10.6 GB |
| B | 128 x 64^2 | 96 x 128.16.16 | 3,145,728 | 1.215e7 / 1.287e7 | 3.98 | 17.1 GB |
| C | 128 x 64^2 | 24 x 128.32.32 | 3,145,728 | 1.359e7 / 1.355e7 | 4.31 | 16.0 GB |
| E | 192 x 32^2 | 24 x 192.16.16 | 1,179,648 | 8.90e6 / 7.95e6 | 7.14 | 9.5 GB |

Readings:
- **Cold start vs developed.** Arm A cold (1.046e7) is within 3.3 % of the developed prod4
  anchor (1.012e7). Below, every arm is scaled by that factor, 0.968.
- **Per-GPU throughput rises with load** up to at least 1.6M zones per GPU: 1.35x from A to C.
- **Blocks per GPU:** at fixed zones 12 large blocks per GPU (C) beat 48 small (B) by 8.5 %: use
  2 x 2 blocks per panel (24 blocks) at every horizontal resolution; re-scan `lb_nmb_eachrank`
  for each new layout (the split is not additive, `tests_rank_imbalance`).
- **Radial cells cost more than proportionally.** E is 20 % lower in zone-cycles/s than A, so
  1.5x radial cells cost 1.86x per cycle. The likely cause is the longer ck chains and the
  O(N^2) twilight beam integral (`tests_rank_imbalance` (a)). It is also the noisiest arm
  (11 % between runs).
- **Memory** ~3.5 GB + 8 kB/zone per GPU (fit to C, D); ~110 GB per APU allows ~12M zones/GPU.
  The 2^14 blocks-per-rank MPI-tag cap is irrelevant (at most 96 blocks here).

**Where the cycle goes.** In the prod4 form, RT is 81.8 % of GPU kernel time and the column
sweep alone is 72 % (`tests_tm_prof_growth/README.md` (7)).
- **Hydro vs MHD.** The MHD-only kernels (hlld fluxes, FC pack, MHD raise-vel, CT) are at
  most ~9 % of kernel time in that profile. hyd4 should therefore be at most ~5-9 % cheaper
  per cycle, at the same CFL dt (section 2). **Not measured.** The hyd4 smoke job 11943303
  (not mine) will give the number.
- **prod3 is not the forecast.** prod3's own links ran at 1.29e7-2.00e7 zone-cycles/s
  (`bench/cs_mhd_prod3/log.out.116112{23,29,30}`). That is 16-25 cycles/s with the older
  plane-parallel binary, 283 rot in 57 h of loop time. The prod4 form costs 1.3-1.9x more per
  cycle.

## 4. Wall time for candidate grids

Formulas:
- **Cycles to rotation R.** Integrated from the prod3 dt history: 8.91e5 cycles to rot 50,
  2.41e6 to rot 150 and 4.43e6 to rot 283. The last matches prod3's restart ncycle 4,429,831.
- **Radial 1.5x.** The cycle count is multiplied by 1.5 for hydro (dt ratio measured at
  t = 0), and by 1.5 x 13.4/11.5 for MHD at max_eta 5e12 (Ohmic, section 2).
- **Horizontal refinement.** The cycle count is unchanged (dt is radially limited).
- **Rate.** cycles/s = arm mean x 0.968.
- **Layout.** All on 1 node / 2 GPUs, 12/12 split, and 24 h apu links.

| candidate | cycles/s | h / rot (steady) | 50 rot | 150 rot | 283 rot | 24 h links to 283 | GPU-h to 283 |
|---|---|---|---|---|---|---|---|
| **current** 128 x 32^2 (prod4) | 12.87 | 0.33 | 19 h | 52 h | **96 h** | 4 | 191 |
| current, 10/14 split | 14.73 | 0.29 | 17 h | 46 h | 84 h | 4 | 167 |
| **1.5x horizontal** 128 x 48^2 | 7.69 | 0.55 | 32 h | 87 h | **160 h** | 7 | 320 |
| **2x horizontal** 128 x 64^2 (32^2 blocks) | 4.18 | 1.01 | 59 h | 160 h | **294 h** | 13 | 589 |
| **1.5x radial** 192 x 32^2, hydro | 6.91 | 0.91 | 54 h | 145 h | **267 h** | 12 | 534 |
| 1.5x radial, MHD (Ohmic 11.5 s) | 6.91 | 1.07 | 63 h | 169 h | 311 h | 13 | 622 |

- **Radial refit at nx1 = 128 (estimate, not measured).** Refit the map to the developed H_p
  profile of section 1.1. That spreads the 128 cells as about 5.5 per H_p everywhere:
  - zone cost unchanged;
  - H_p-minimum cells 2.9 -> ~5.5;
  - dr at the binding radius 528 -> ~265 km, so dt ~10 s and about 2x the current wall time
    (~190 h to rot 283).
  - That is cheaper than nx1 = 192 for a larger gain at the H_p minimum. The coarser top and
    base cost nothing that the requirement protects.
- **2 nodes / 4 GPUs.** Only worth it for the 2x-horizontal grid, where each GPU would hold
  786k zones. The off-node efficiency at that load is unmeasured. If it matches the 08-26
  1.92M-cell point (2 -> 4 GPUs gave 1.52x), 294 h becomes ~195 h on 2 nodes. At the current
  grid 4 GPUs gave only 1.28x (`tests_mpi_gate_0922` g2a -> g4), so stay on 1 node.

**Physically sufficient duration** (prod3 `dhj.mhd.hst`, mean in +-2.5 rot windows divided by
the rot 183-283 mean):
- total KE: 0.91 (rot 10), 1.16 (20), 1.26 (40), 1.16 (100), 0.95 (150), 1.05 (200);
- radial KE: 0.56 (10), 0.82 (40), 0.97 (150);
- ME: 0.55 (10), 0.58 (40), 1.28 (100), 1.07 (200), 0.90 (250).

The upper jet saturates by rot ~40 on sp, but the deep equatorial westward flow strengthens
monotonically to rot 164 (memory `dhj-jet-shallow-westward-deep`).
- **~50 rotations** suffice for the upper-atmosphere observables (jet, day-night contrast,
  terminator T(p)).
- **At least 150-200** are needed for the deep flow and ME.

## 5. Recommendation

1. **prod4 as launched is sound for cost:** 96 h to rot 283 on 1 node, 4 links; 84 h with
   the 10/14 split.
   It delivers 5.0 median / 2.9 minimum cells per H_p (not 10) and 2-4 cells across the jet.
2. **Next horizontal grid: 48^2 per panel**, blocks 128 x 24 x 24 (24 blocks, 12 per GPU),
   1 node.
   1.67x per rotation, dt unchanged: 160 h to rot 283 (7 links), 32 h to rot 50. 64^2 costs 3.1x
   (294 h, 13 links); only on 2 nodes, after measuring the off-node efficiency at that load.
3. **Radial:** refit the stretch to the developed H_p profile at nx1 = 128 before adding cells.
   The requirement "10 per H_p everywhere below 1e-6 bar" needs about 191 cells with a matched
   map, i.e. nx1 ~ 210-220 with a coarse top. With dt ~ dr/(|v_r|+c_s) at the H_p minimum,
   that costs about 3.6x per rotation on dt alone (dr 528 -> 145 km at the binding radius)
   plus about 1.9x per cycle (estimate, extrapolated from arm E).
   - The requirement's physical necessity is **untested** in a developed state.
4. **For MHD with any radial refinement**, scale max_eta with dr_min^2 or enable use_rkg_sts.
   Otherwise the Ohmic cap binds (11.5 s at nx1 = 192).

**The one test that would reduce the uncertainty most:** a developed-state radial resolution
pair.
- Arms: the hyd4 twin from scratch at nx1 = 128 vs nx1 = 192 (same map), both to rot 20;
  about 10 h + 27 h on 1 node each (section 4 rates; 4.42e5 cycles to rot 20).
- Compare at matched rotations: terminator T(p) and 1e-6 bar radius, H2-front pressure, jet
  u_eq and FWHM at 1e-4..1e-2 bar, total and radial KE.
- Reading: a difference beyond the run-to-run spread (seed ensemble, memory
  `cs-vertex-dt-collapse-0907-defaults`) says the radial grid, which today delivers 2.9 cells
  per H_p at the H_p minimum, sets the answer.
- A cheaper first check needs no GPU: redo the jet-width table above on a common setup, i.e.
  the sp arm with the cs RT switches suggested in memory `cs-vs-sp-dhj-rot50`.
