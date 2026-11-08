# MeshBlock-layout invariance suite

Checks that results do not depend on the MeshBlock decomposition or on the MPI rank count.
It was written after the cubed-sphere seam-halo bug of 11-03, in which the seam stencils were
clamped to the sending MeshBlock instead of the panel. That bug went unnoticed for months
because conservation still held.

## How it works

Each case in `cases.py` runs the same problem several times from identical fresh analytic ICs:

- the **reference layout** (the first one, normally one MeshBlock per panel);
- the reference layout again with a **1-ulp perturbation** (default `time/cfl_number` moved
  to the next double);
- every **other layout**: more blocks per panel, a phi or theta split, a radial split,
  1 to 8 MPI ranks.

Every run writes double-precision `.bin` output (`bin_double`, `dcycle = 1`).
`libin.py` assembles the MeshBlocks into global `(panel, k, j, i)` arrays. Cubed-sphere gids
are panel-major. The comparison is done cell by cell for every listed variable at cycle 1 and
at cycle `ncyc`.

The 1-ulp run measures how far plain round-off spreads by that cycle. The verdict for each
(layout, cycle, variable) is one of:

| verdict  | relative max difference vs the reference layout            |
|----------|-------------------------------------------------------------|
| BITWISE  | 0                                                           |
| ROUNDOFF | <= max(10 x ulp-ref, 100 eps)                               |
| SUSPECT  | in between                                                  |
| FAIL     | > max(1000 x ulp-ref, 1e-10)                                |

For SUSPECT and FAIL the report also gives the cell with the largest difference, and the
fraction of differing cells that lie within `ng` of an INTERNAL block end of the test layout.
A value near 100 % is the signature of this bug class. For scale, the 11-03 clamp gives
about 1e-5 after one cycle, where round-off is about 1e-16.

A case can carry `expect = 'FAIL'`. This keeps a known layout dependence visible: the
`cs_seam_blockend = clamp` default and the r1m clamp deck are both kept this way. The summary
flags any case whose outcome differs from its expectation.

## Running

```
# bins.cfg: "family = /path/to/athena" for bip (built_in_pgens), cs (cs_test), sp (sp_test),
#           dhj (deep_hot_jupiter_rt), he (he_star_m1), acc (ry_per_accretor); CPU builds, MPI on.
python3 li.py --list
python3 li.py --bins bins.cfg --work /scratch/li_runs [--only case1,case2] [--np-max 16]
python3 li.py --work /scratch/li_runs --compare-only     # re-compare existing runs
```

Exit status is 0 when every case meets its expectation. The report is printed and written to
`<work>/REPORT.txt`.

On Caltech, use `run_li.sh`. It loads gcc 13.2 and openmpi 5.0.1. Do not set
`OMPI_MCA_io=romio321`: openmpi 5 has no such component. The binaries are listed in
`bins.caltech.cfg`. Launch ON a compute node of a held allocation, detached:

```
ssh <node> "setsid nohup bash -l .../run_li.sh <work> --only ... > <work>.out 2>&1 < /dev/null &"
```

The dhj cases need the WASP-121b packages, the ck tables and the IC profile (env
`LI_DHJ_PKG`, default `/resnick/groups/carnegie_poc/jingze/w121r1_1101/pkg`). The whole suite
takes about 15 min on 16 cores; the cubed-sphere and Cartesian cases take seconds each.

## Coverage

- **Cartesian** (`built_in_pgens`): 3-D hydro and MHD linear waves (CT edges shared by 4 or 8
  blocks), 2-D implode with reflecting walls (physical-BC corners), Orszag-Tang, and an M1
  beam with user/vacuum boundaries.
- **Cubed sphere** (`cs_test`):
  - hydro and MHD blasts through a seam and a cube vertex, with clamp and with extend;
  - nghost 3 with ppm4, the HG/RSG halo depth;
  - the r1m seam keys (`cs_seam_flux = shared`, `cs_seam_mhd = shared_o`, `uct_hlld`);
  - Ohmic plus ambipolar diffusion with RKG STS;
  - scalar diffusion and a passive scalar;
  - probes at C8/C16/C32 that separate block size from blocks per panel.
- **Spherical polar** (`sp_test`):
  - hydro and MHD blasts across the polar boundary, with the phi direction split into 2, 4
    and 8 blocks, a theta split, and a radial split for hydro;
  - `polar_ring_avg` (auto group sizes), the ring Fourier filter, `polar_emf_diss` and
    `polar_hlle_rows`;
  - resistive decay with Ohmic, ambipolar and RKG.
- **dhj** (`deep_hot_jupiter_rt`), using the production r1m v6 MHD deck and the hydro planet
  deck on C8 x 32:
  - physics: well-balance cache, top sponge, Roche tide plus the `wb_tan` probe, `outflow1`
    and `hse` top BCs, bottom sponge, floors;
  - ck: the implicit ck two-stream (spherical beam) at cadence 4 and 1, plus a
    `ck_impl_jreuse = 0` control;
  - MHD: RKG Ohmic plus ambipolar diffusion and tracers;
  - layouts: 1 vs 2x2 blocks per panel and 1, 3, 4 or 6 ranks. dhj refuses a radial split.

- **HG** (he_star_m1 production deck, C16 x 704, M1 with vet_col, vet_gd twin and vet_sc),
  comparing 1 vs 2x2 blocks per panel and 3 vs 6 ranks. The 1-ulp seed is `he_gm`, because
  dt is set by the `vet_gd_rotate_dt` alignment, not by the cfl. This deck is so sensitive
  that 1 ulp already gives 1.5e-6 in rho after one cycle, so the case resolves layout bugs
  only above that level.
- **sp accretor** (ry_per_accretor production deck, 500 x 48 x 64): phi split into 2 or
  4 blocks, with vet_gd on and off.
- **red_giant**: the vpert seed with and without `vpert_global`.

- **merge_lat** (`<hydro>/merge_lat_blockend`): f 2, 4 and 8 and the RSG pattern 8:4:2 with
  the RSG E key set (nghost 3, fofc, `cs_seam_flux = shared`). Each runs as flat (the old
  rule, expect FAIL) and as exchange. The same comparison is done on the RSG E production
  deck at C32 (`rsg_mlat_x`).
- **sp polar rows**: ideal UCT and GS07, Ohmic and ambipolar, each with a 1-ulp IC
  reference (`problem/b0`).
- **Cube-vertex probe** `r1m_v8_vertex_cc`: the r1m v8 deck with the vertexacc-1103
  `cs_vertex_cc_perturb` diagnostic (bin family `dhjvtx`).

## Not covered yet

- Restart starts from a permutation remap, i.e. production states. Fresh analytic ICs only.
- Refined meshes: `libin.assemble` handles single-level meshes only.
- GPU builds: run the same cases with GPU binaries in `bins.cfg`. The suite is
  backend-agnostic, but it uses mpirun on one node.
- The ulp reference may not cover the real round-off seed:
  - Cases whose dt is not set by `time/cfl_number` need their own `perturb` key, otherwise
    the ulp reference reads 0. Check the `ulp-ref` line.
  - Cases whose IC is not bitwise identical across layouts need an IC-level perturbation.
    An example is a field built from a vector potential on per-block coordinates, which
    differs in the last bit. A cfl perturbation then underestimates the round-off seed,
    and the case reads SUSPECT, or FAIL in GS07 near v = 0.

## layout_invariant and the two gate kinds (layinv-1105, USER 11-05)

`<mesh>/layout_invariant` (fresh-start default true; a restart whose input lacks the key keeps
false, the old code paths, bitwise) makes a run independent of the MeshBlock layout and of the
rank count:

- every global or shell SUM that feeds the state uses the exact, order-independent accumulator
  `src/utils/repro_sum.hpp`: the M1 BiCGStab dot products, the vet_gd twin / half-range /
  trans shell means, the ic_rad shell sums and L_top, the he_star IC shell sums
  (he_ic_rad_atm_post, he_ic_equil) and the adaptive MLT closure;
- the spherical-polar theta and phi faces/centres come from the GLOBAL cell index
  (`<mesh>/sp_theta_global`, `sp_phi_global` also on wedges); a stretched theta is applied to
  the global uniform coordinate exactly as before;
- the vet_gd sweep takes its (theta, phi), spacings and fractional upwind indices at the
  global magnitude (`src/rad_m1/rad_m1_vgdlg.hpp`).

The M1 PRECONDITIONER can still depend on the layout: `rbgs`, `rbgs_fwd`, `mg`, `mg_gc` and
`mg_gf` are block-local (red-black parity and couplings inside the MeshBlock; mg coarsens inside
the block).  Only `implicit_precond = line` (x1 lines, blocks span x1) is layout-independent.
With a block-local preconditioner the Krylov / Picard iterations stop at their tolerances on
different iterates, so two layouts agree to the TOLERANCE, not to round-off.  Hence two gates:

| gate      | decks                                                     | case keys                  | pass                     |
|-----------|-----------------------------------------------------------|----------------------------|--------------------------|
| BITWISE   | line-preconditioned decks, hydro-only cases, the cycle-0 IC | `require_bitwise = True`   | every variable identical |
| TOLERANCE | mg / rbgs decks (production R136a1, BSG, accretor, HG, RSG) | `tolref = {key: factor}`   | layout difference <= `tol_factor` (default 1) x the spread of the reference layout re-run with the tolerances scaled by the factors (e.g. `implicit_lin_tol`, `implicit_lin_cnorm` (the linear stop test when named), `implicit_tol` x 0.5) |

`tolref` replaces the 1-ulp perturbation of the reference run (`<layout>_ulp`) by the scaled
tolerances; the verdict is BITWISE, TOLERANCE or FAIL.  USER 11-05 chose this split (option a): no
preconditioner change.

**Gate (a1), the one in use (USER 11-05, after the tolerance reference proved too tight).**
The linear solves overshoot their stop tests, so the x0.5 tolerance re-run barely moves (dens
2e-19, F 4e-13 on the BSG deck) while a different block-local preconditioner lands on a different
iterate (dens 4e-14, F 2.5e-12, the same as the old binary's layout pair).  A mg / rbgs deck
therefore PASSES when (1) its layout pair is no worse (x2) than the OLD binary's same layout pair
(`ref_old = True`; bins family `<family>_old`; verdict NOWORSE) AND (2) the same deck with
`implicit_precond = line` is BITWISE across the layouts (`require_bitwise`, e.g. `acc_sp_th_line`).
The `tolref` machinery stays for reference.

## Restarts of pre-1105 decks: name `mesh/layout_invariant=true` (layoutdet-1107)

A restart whose input lacks `<mesh>/layout_invariant` keeps `false` (old block-ordered sums).
The HG and RSG production decks are such restarts. Their rank-count dependence (1 vs 2 vs 4
ranks, same 24 MeshBlocks: dens ~1e-13, eint ~1e-11, m1_e ~1e-10 after 4 cycles; near-zero
velocities flip sign) comes from these old sums. The layout det trace below puts the first
difference at the first solve of the cycle in `iw(M1_IW_EN)` (E^n + dt esrc) while u0, opac and
the fluid state are still identical: `esrc` comes from the adaptive MLT closure
(`HsAdaptiveUpdate`, he_star_m1.cpp), whose radial shell sums are a serial loop over the
rank's own MeshBlocks followed by the rank-ordered `HsDetSum`, so the association of the terms
changes with the rank count. The BiCGStab dot products (`ImplicitBiCGStabRun`: a Kokkos
reduction over the rank's pack, then `M1GlobalSum2`) have the same defect.
Adding `mesh/layout_invariant=true`
to the restart command line (it is then recorded in the next restart) makes them bitwise across
rank counts, with no code change (prod-1107v 698b82413, layoutdet-1107 gates of 11-08):

- hgb and rsb small C32 states, control and r8-cs + T: 2 vs 1 rank (MI300A), 4 ranks on 2 nodes vs
  1 rank (MI300A), 4 vs 3 ranks (A100): bins and restart data identical;
- restart (R+2 -> rst -> +2) identical to R+4, and run-to-run identical;
- iteration counts unchanged; cost +8 % per cycle on hgb (2 MI300A ranks, 20 cycles, interleaved).

This covers a change of the RANK COUNT with the same MeshBlocks. A different MeshBlock layout
still differs at the Krylov tolerance with the block-local preconditioners (`mg`, `rbgs`, see
above).

Diagnostic: env `M1_DTRACE_LAYOUT=1` with `<rad_m1>/implicit_det_trace = N` prints on rank 0
layout-comparable hashes (global position mixed, exact uint64 sums over all ranks) of the
solve's arrays at clos_in, opac_in, isolve_in, solve_in, pass_pre, vet_built, pre_bcg, post_bcg and c2p_in for the
first N solves; diff the `<dtl>` lines of two layouts to find the first array that differs.
