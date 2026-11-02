# NOTE 2026-11-01 ~20:40 PST Caltech -> viper (USER): fused Strang + ck Strang results, order ladders running; CPU non-determinism found

Seen: 69f87070 (fin51 = 7f0e3dd9, halofc merged, thanks), c586f965 (ck_impl_maxit). USER: KEEP ck_impl_maxit = 16 (deck unchanged).

## 1. Branches (pushed, not merged)
- rkgfuse-1102 @ 45c4bf29 (tests ran on 8e7e35eb; based on hfcstrang fcc40f98 = halofc + rkgstrang):
  - <mhd>/rkg_strang_fuse: the closing RKG half R(dt_n/2) is fused with the next opening half into ONE super-step over
    (dt_n+dt_{n+1})/2 at the start of cycle n+1, after before_timeintegrator. Halves are unfused (synced) before any
    output/rst, at loop end, and with AMR (fused FATALs with AMR). The pending half is restored on dt_retry.
  - <problem>/ck_strang: the ck RT runs as two dt/2 halves nested inside the RKG halves (R ck H ck R). The cadence counts
    cycles: only the first half of an ncycle % ck_impl_every == 0 cycle is a full solve.
  - <mhd>/tracer_strang (45c4bf29, gates running): TracerSource in two dt/2 halves (R Tr ck H ck Tr R).
  - All keys default false and are read only when named.
- rkgfuse-fin51-1102 @ b1954614: port onto 3957c956 (ck_strang edits re-applied in the split TUs). It will be moved to
  7f0e3dd9 after the order ladders.

## 2. Gates (H200): all IDENTICAL
- key absent, 40 cycles, rkg_strang off and on vs fcc40f98: rst data IDENTICAL
- 26+26 vs 52 (dump forced at 8800) for fuse / ck_strang / full (fuse + ck_strang): IDENTICAL

## 3. Cost, r1m, 1 H200 (halofc included), interleaved x2

| arm | ms/cycle (diag window) | RKG stages/cycle | kernel ms/cycle |
|---|---|---|---|
| lie | 74.4 | 21 | 66.9 |
| strang | 78.0 | 30 | 84.2 |
| fused | 62.2 | 21 | 67.1 |
| full | 64.6 | 21 | |

- Fused does the same GPU kernel work as lie. Its wall lead is 7-11 % by whole-run zone-cycles/s, and 16 % only in the diag
  window, which caught one slow lie interval. The cause is host-side (sync/overlap placement), not yet explained.
- Hydro W121 p1, ck_strang: 27.18 -> 28.94 ms/cycle (+6.5 %). Of that, ~0.2 ms is the extra linearised half and the rest is the
  RTOpSplitBvals tail (halo + BCs + c2p) after the extra split call. rt_pre_opac/pre_tp/apply run only in full calls; Newton
  passes per full call drop 3.9 -> 3.0 (dt/2).
- New here: an energy-only post-split refresh (branch ckbval-1102 off 7f0e3dd9, bitwise only, in progress). Expect -4-5 %
  on hydro dhj, about 2x that with ck_strang.

## 4. Accuracy (0.02 rot from off_0p50, your c8017e87 spec, strang4 = reference)
- fused vs strang4: equal to strang within 2 % above 1 bar. At p > 1 bar the B error is ~2.7x strang (lie level, ~1e-5 of the
  field).
- full ladder: (full - full16)/(full4 - full16) = 1.2-3.0 by band and field (2nd order = 17, 1st = 5). At 1e-3..1 bar and
  p > 1 bar these differences are 10-5000x the seed noise, so the poor convergence is real, not chaos. The remaining error
  does not shrink with dt. Prime suspect: the cycle-counted cadences (ck_impl_every 4, ck_impl_xstep 8, ck_vef_every 150 =
  430 ms/call). The error attribution ladders are running:
  - (a) all cadences 1;
  - (a') ck cadences 1 only, VEF 150;
  - (b) (a) + tracer_strang;
  - each at cfl / cfl/4 / cfl/16, plus a seed pair at cfl/16.
- Counters: no eos_fail / vceil / fofc / NaN in any arm; floors per c2p call not above lie; ck cap/guard/fallback 0.
  One flag: hydro eos_tclamp 232 -> 380 under ck_strang (x1.64 vs x1.33 more c2p calls).

## 5. CPU non-determinism (please note for your CPU gates)
A single-rank Serial dhj-tree CPU build (PROBLEM = a local LinearWave wrapper) is NOT reproducible run to run: ~1e-8 absolute
from cycle 3, independent of amplitude, also in ideal MHD and with Lie. Our accretor sp40 CPU runs were reproducible only
with OMPI_MCA_coll_hcoll_enable=0 OMPI_MCA_coll_ucc_enable=0. A diagnosis agent is running: HPC-X collectives, OpenMP,
stock linear_wave vs the wrapper, valgrind/initcheck. Your clean CPU convergence (Lie 1st, Strang 2nd) stands; ours was
noise-limited, so we rely on yours plus the GPU gates.

## 6. Other
- HG global switched to vet_gd twin OFF at 20:18 (A/B passed); RSG has been on twin off since 17:09.
- cs accretor: <rad_m1>/time2_floor_eint = 1e-7 (csacc-1031 20b3b78d) cures the BE fallbacks (465 -> 2 BE steps, 1.64x
  faster, L/Min/Mnzf within ~1e-4). It is being adopted in the decks; accwbtan-1101 is being merged into csacc-1031. A
  csacc -> fin51 test-merge branch + NOTE to you will follow (you own the fin51-r-fix merge).
- r1m: NOT launched. The user decides the scheme after the order ladders.
