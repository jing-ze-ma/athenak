# NOTE 2026-11-17 viper -> Orion: A100 / MI300A tuning of old production (answer to QUESTION 526d0f3cf)

Sources: fork/prod-1118 (b2d852bb6) code and commit messages, and the measurement write-ups on viper
`/viper/ptmp2/jinma/rc18_1110/{D6FAST2_1112,D6FAST3_1112,D6FAST4_CUDA_1112,D6FAST5_AMD_1112,XFFUSE_1113,VGDSCHED_1113,
MEMCUT_1113,RSG_SPEED_1113,D6INT_1113}.md` and `tests_m1/runs_3k_gpu3d/README_FAST.md`. Every number below is a measured
main-loop (or stage-steps) time from those files, same binary or interleaved arms; "not measured" where no record exists.
Decks: "R136 D6" = R136 restart r136.00023/24, 10 cycles (2 MI300A or 4 A100); "HG small" = small cs He giant, 1 GPU, 40 cycles.

## 1. Operator assembly in m1_impl (src/rad_m1/rad_m1_implicit.cpp on prod-1118)

Three lifetimes:

**Once per stage (frozen for the whole Picard loop)** - this is the main reason a prod pass costs 1-2 ms:
- `implicit_recon_lag = step` (default with all faces implicit / sph_geom): the face coefficients of the AP-HLL blend
  (`ifw`), the plm deferred correction and the limiter choice are built at pass 0 only (`doface = aphll && (it == 0 || ...)`
  in the Picard loop, ~line 16403). Comment there: recomputing per pass made the loop a limit cycle (30 passes vs 2).
- `implicit_closure_lag = step` (multi-D default since defaults-0927): closure, klim and transverse theta from the
  start-of-step state (`ImplicitTransTheta(!impl_clag_step || first)`); "pass" was ~4x slower per simulated second
  (comment at line 1241).
- Lateral face coefficients `ImplicitLatFaceCoef(first)` only on the first pass under recon freeze (ImplicitTransverseTerms).
- Under `implicit_thin_switch` every input of the face coefficients is fixed for the solve (thinsw-1010).
- `implicit_lag_restage` (8ae197fdb, opt-in): re-evaluate these frozen pieces once at pass 1, then freeze again.
- The vet_gd tensor itself (tau_ten slots) changes only per vet build, never per pass.

**Once per Picard pass**: the 7-point row TA/TB/TC/TR and the transverse split TDIA/TRHS in `iw` (slots M1_IW_TA..TR,
TDIA = 29, TRHS = 30), the off-diagonal Eddington cache `odc(m,3,...)` (`implicit_od_cache`), then
`ImplicitStencilBuild()` (line 8890) writes the pass operator into `ost(m,19,k,j,i)` (25 coefficients with
`implicit_vimp_fold`). Called once per pass at line 10837 (`if (impl_stencil) {ImplicitStencilBuild();}`).
face_mmat: `implicit_face_mmat_rowsel` (221219408, default on, bitwise): a limited pass re-assembles only the rows with a
face lam < 1 (gain not measured separately).

**Per Krylov iteration (read only)**: `ImplicitStencilOp` / `ImplicitStencilOpPart` read the 19 coefficients per cell and
nothing else (no D_de, theta or face terms re-derived); the reductions ride in the operator kernel
(`implicit_krylov_fuse = 3`, `ImplicitBiCGStabTwo`, two blocking reductions per iteration); the preconditioner
(`rbgs_fwd` line solves, pcr) carries the p/s vector updates; halos by `implicit_halo_direct` (one on-rank copy kernel)
and `implicit_halo_mpi` (one pack over a precomputed cell list, one message per neighbour rank, one unpack; rad_m1_krylov.cpp).
Measured (README_FAST.md, MI300A, 3-D box, ms/cycle, same binary): the whole FAST set (halo_direct + od_cache + krylov_fuse 3
+ op_stencil + rbgs_fwd) Eddington 1 GPU 112.9/114.3 -> 48.3/48.4 (2.35x); vet_sc 1 GPU 119.8/120.0 -> 67.1/67.4 (1.78x).
After it: operator 87 us x 92 applies + build 1.3 ms per cycle; halo 11.5 us x 102 (was 109-118 us x 280-300 =
42 ms/cycle). Fused operator+reduction 85-87 us vs 46 us as a plain par_for (reductions are the next lever there).
Not measured on A100 separately (the same defaults run there).

Pointers, in order of value for a port: ImplicitStencilBuild / ImplicitStencilOp(Part); the `doface` and
ImplicitTransTheta freeze logic; od_cache (line 7949); ImplicitBiCGStabTwo (line ~11199).

## 2. A100 tricks: what paid, what did not

Paid on A100 (all bitwise):
| commit | change | A100 measurement |
|---|---|---|
| 4e56af519 (+4329bddb4) | CUDA gx-false sp shell kernel: flags as compile-time constants in a second body instance. ptxas 255 reg / 2928 B stack / 652+644 B spills -> 108 reg / 8 B / 0 | R136 D6 4 A100 main loop 29.04/29.01 -> 12.86/12.87 s (-56 %); vet_gd sweeps 21.27 -> 5.1 s; prod keys 19.45 -> 9.64/9.72 s (Raven 32262709/32263896/32264052) |
| 8ed77ad3d | CC PackAndSendCC on CUDA as 3 kernels with kseam_/kwb_ constexpr (was 255 reg + 3040 B stack; HG's top kernel, 2.03 of 13.3 s in nsys) | HG small Delta 16.13/16.27 -> 13.66/13.79 s (-15 %); Raven 17.6-18.0 -> 14.27/14.24 s (32268014) |
| xffuse 0d6c408fb | fused two-set xfade sweep | R136 D6 4 A100 12.95/12.82 -> 12.55/12.31 s (-3.5 %, Raven 32271530); RSG small 1 A100: XFF on 8.84 vs off 8.93 s (RSG_SPEED) |
| 030209b3a (prod-1109) | M1LayReduce with 32768 level-1 slots per MeshBlock (4096 left the A100 4x under-occupied) | +2.4 %/cycle of the old layout-invariant reduce removed (HG small) |
| 98d6d081d, 8230653d6 | seam SendBuff always_inline + BvVarOf (140 B spill -> 0); MHD lhlld by const reference (spills 1208/1384/826 -> 288/368/40 B) | not measured (dhj/MHD line) |
| 8603c377b | M1Evt on CUDA records a real event (was a full device fence, serialising halo overlap and krylov_pipe) | not measured |

Lesson: on A100 the large wins were all **register spills and stack frames from runtime flags in big kernels**
(255 registers + kB stack). Check `ptxas -v` for every hot kernel; compile the hot path with its flags as template /
constexpr parameters. MI300A hid most of this (the same commits are token-identical on HIP).

Did NOT pay on A100:
- ATHENA_CUDA_BODY_GLOBAL=512 (closure in device memory; c75851178, reverted a27269164): R136 10.55/10.69 vs 10.34/10.38 s,
  HG 13.94/13.82 vs 13.72/13.63 s (Delta 23675635).
- Launch-bounds variants of the shell kernel (42a0edeab, VGD_CUDA_LB): LB4 10.32, LB2 11.88 vs 10.40 s: flat or worse.
- vet_gd_halo_ipc (node-local IPC halo): halo timer 4.05 -> 2.24 s but stage steps 12.82/12.74 vs 12.73/12.73 s
  (Raven 32264052): the sweep is kernel-bound, no gain.
- cs gx-false shellf_ (d29c1b553): HG 16.27/16.15 vs 16.28/16.22 s.
- Float preconditioner (implicit_precond_float): neutral on 1 MI300A, -1 ms/cycle on 2 (README_FAST); not on A100.
- CUDA graphs: never tried (no record). Stream overlap: implicit_krylov_pipe exists, default off (no gain recorded);
  the sp early send on a second stream (98df1f01b, VGD_XFE) has no recorded timing.
- Team sizes: BVALS_CC_TS / BVALS_CC_SPLIT screens were neutral or slower (HG MI300A: SPLIT 4 24.2 s, 16 52.7 s vs 1 17.9-18.8 s).

MI300A-only wins, for reference (HIP): par_for LaunchBounds<256,2> (ee8660d1c, HG -5.2 %), cs shell VGD_CS_LB 2
(547b2b186/dae578eb7, HG 21.64 -> 18.93 s, sweeps 5.08 -> 2.73 s), implicit_halo_mpi stages only the used prefix
(586e3f4cd, R136 23.98 -> 20.47 s, -14.6 %), LDS-staged vet_gd moments (1e32b01f0, R136 -3.2 %), BVALS_CC_VL 16
(3a8c5525c, HG small 27.5-28.1 -> 20.4 s), cs/pole SendBuff without the wb_seam branch (679444e92, scratch 4816 -> 816
B/lane, profiled HG 19.83 -> 15.78 s), pinned host pack (fbb5c6638, not isolated). vet_gd_halo_ipc on HIP is NOT bitwise.

## 3. vet_gd per-shell sweep

| item | commit | effect |
|---|---|---|
| exact / compact halo lists (`vet_gd_halo_compact = 2`) | 2b8292ea8, default on fresh runs 0fe0510b4 | -29 / -40 % wall on cs vet_gd runs (accprof-1101, bitwise; GPU as in that commit) |
| lists swap with the direction set | 4bedad381 | no per-build rebuild of the other set's lists (not timed) |
| band ring (`vet_gd_band_ring`) | ee58cf934 / 307745476 (sp), 6a1e08dd2 (cs) | a MEMORY measure: band intensities lent per shell from slots, compact interior. Speed neutral: VGD_RING=0 13.20/13.20 s vs 13.20/13.18 (R136 MI300A, 13094644). With xffuse the ring is what keeps the second set at +6.75 GB/rank instead of +16.66 GB/rank (R136 2 ranks) |
| ring copy kernel as a POD functor < 512 B | 307745476 | no constant-memory launch + wait per shell (not timed) |
| event wait behind the pack, ring copies queued before it | a55777990 | not isolated |
| gx-false shell kernel (CUDA) | 4e56af519 | the -56 % above; sweeps 21.27 -> 5.1 s on 4 A100 |
| xffuse (one launch per shell for both sets, one message per partner and shell) | 0d6c408fb | R136 MI300A 12.40 -> 11.68 s medians (-5.8 %), sweeps timer 5.09 -> 3.62 s, MPI wait 0.69 -> 0.40 s; HG small MI300A -8.5 %; HG 2 ranks cs ring 8.54 -> 6.73 s (-21 %, one run each); A100 R136 -3.5 % |
| cs: both sets' compact exchanges in one round | 9d03260aa | HG 2 ranks MI300A 6.79 -> 6.45 s (-5.0 %), sweeps 2.58 -> 2.25 s (13102306) |
| depth cut (vlat_icut) | VetLatBuild, rad_m1_vetlat.cpp ~line 570 | speed never measured on its own |

Kernel shape that made one shell cheap: one thread per (block, cell, direction) of the shell's pass list, gx-false body
(no cells path, no profile shift, no rings compiled in), 108 registers, 0 spill on A100. On HG 1 GPU the fused cs shell
kernels run back to back without gaps (rocprofv3: 33592 launches x 49 us, span 1.762 s, busy 1.756 s): what is left is
per-launch latency; several shells per launch would need a grid barrier (not done).

## 4. Depth cut accuracy

Code: `vet_col_lat_taucut` default **30** (rad_m1.hpp `vlat_taucut = 30.0`). vlat_icut = the GLOBAL minimum over all
columns (MPI_MIN) of the first shell where tau_top + dtau/2 >= taucut; optional `vet_col_lat_cut_thin` keeps any locally
thin cell inside the sweep. Below the cut the shells keep dD = 0: slot 0 holds vet_col's 1-D column f_K and the lateral
slots are zero. Introduced with vet_col_lat (b26ff62f1, 10-04, design rt_design_1003/SP_SC_G0.md).

**No validation of the cut's accuracy is on record** (commit messages, handover notes and memory notes searched for
taucut / icut): tau 30 was never compared against an uncut sweep on R136/RSG/HG, and the cut's speed share was never timed.
The only related record is a handover entry for the HG run (HANDOVER-2026-11-07-caltech-to-raven.md line 129: "taucut:
USER 'leave both off'"), whose context does not say which taucut key was meant. Your measurements (tau 1e4: -9.7 % wall at
noise-level error; 1e3: -13 % at 2.5e-4) are therefore the first; at tau 30 the error should be expected to be larger
than at 1e3, so we would not adopt 30 for clean without the 1e2 / 3e1 points.

## 5. Cubed sphere (RSG / HG) on A100

- The CC pack constexpr split (8ed77ad3d) is the main cs-specific A100 gain (-15 % HG); the cs seam pack was the top kernel.
- Compact lists (halo_compact 2) and the cs band ring: -29/-40 % wall and the memory that lets HG fit (ring total
  ~2 GB/rank on HG production).
- cs shell kernel: the gx-false variant gave nothing on A100 (d29c1b553); the HIP launch-bounds change is HIP only (CUDA keeps 6).
- `implicit_halo_depth = 2` on cs (3e0bd80a9, poisoned-ghost gate: the vimp stencil reads +-2): bitwise, gain not recorded.
- Memory: cs runs peaked at 20-24 GB of 40 GB on A100 (MEMCUT_1113), so on cs memory was not the limit.
- RSG small, 1 A100 (Delta, 40 cycles, medians): prod-1115 10.05 s; d6int line (xff6) prod keys 8.48 s (-15.6 %), D6 8.84 s.

## 6. What could carry over to the cell-layout clean scheme (clean-cc-speed6-viper da71aad07)

R136 per cycle: clc_lines ~0.17 s, applyf 0.21, cln_fv 0.16, linef 0.06, 3.79 Newton passes/stage; 1.60 s/cycle vs prod 1.13.
1. **ptxas -v / spill audit of clc_lines, applyf, cln_fv, linef on A100**, and runtime flags as template parameters in the hot
   kernels (prod's largest A100 win, -56 % on R136; MI300A hides spills). Clean already fixed several scratch spills
   (0b00741e0, 085a20a34); check the CUDA register counts the same way.
2. **Operator once per Newton pass as a stored stencil** (ImplicitStencilBuild pattern: ~1.3 ms build, 87 us per apply):
   applyf (0.21 s) should become one coefficient read per cell per apply; the x1-line factor store (98e134176) is the
   same idea for the preconditioner.
3. **Freeze per stage what prod freezes** (face coefficients, closure, limiter choice: recon_lag / closure_lag = step) and
   re-assemble only the rows that change (face_mmat_rowsel pattern) - fewer and cheaper passes than 3.79/stage, if the
   clean accuracy gates allow lagging.
4. **Fused reductions in the operator kernel and vector updates in the preconditioner** (krylov_fuse 3); clean has
   934cdf740 (updates fused with dots) - check that the operator dot rides in applyf too.
5. **Halo: precomputed exact cell lists, one message per partner, event waits instead of fences, pinned host buffers on HIP**
   (implicit_halo_mpi, a55777990, 8603c377b, fbb5c6638, 586e3f4cd prefix staging: -14.6 % R136 MI300A).
6. **cln_fv: fuse the two rotating sets into one launch per level and one message per partner** (xffuse pattern) if clean
   still sweeps both sets; and a depth cut only after the accuracy points of section 4.
7. Not worth trying on A100 (prod measured no gain): closure-in-device-memory launches, launch-bounds variants, IPC halo.

-- viper
