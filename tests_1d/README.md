# he4_presn 1-D gates (2026-09-17, viper, branch `he4-presn-global`)

What this directory is: the artefacts of the `problem/ic_profile` reader gate and of the
G.3/G.4/G.5 1-D column gates for the 4 Msun presupernova He star on the cubed sphere.
All runs are `-DPROBLEM=red_giant`, GPU (HIP, `Kokkos_ARCH_AMD_GFX942_APU`, MPI on),
built into `../build_gpu_rg`.  All `.bin`/`.rst` have been deleted; `.hst`, the logs and
the cycle-0 column are kept.

Build:

    module purge && module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
    cd ../build_gpu_rg && make -j16

The 1-D configuration is the production input with the angular grid collapsed to one
cell per panel and the seed off:

    -i ../inputs/hydro/he4_presn_cs.athinput \
       mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0

## C. Regression: the reader must not change the march path  -- PASS

`sub_reg.sh`: `inputs/hydro/red_giant_cs.athinput time/nlim=100`, run twice on one GPU,
first with the binary built from HEAD (`2b5a1684`, saved as `athena_ref`) and then with
the binary carrying the reader.  The input's three table paths are orion paths and are
overridden on the command line.

    cmp reg_ref/rg.hydro.hst reg_new/rg.hydro.hst   -> BITWISE IDENTICAL
    cmp reg_ref/column.txt   reg_new/column.txt     -> BITWISE IDENTICAL
    (md5s kept in regression_md5.txt; the 1 MB column.txt pair was deleted for quota)

## D(i). Cycle-0 column read-back  -- PASS

`column_he4_ic_cycle0.txt` (the run's own `column_dump`), 4080 fine-grid nodes,
against `bench/hestar_presn/ic_he4_presn_sph.txt` and `column_he4_presn_sph.txt` over the
mesh `1.18585e11 .. 2.4057e11 cm`:

| quantity | max relative difference | where |
| --- | --- | --- |
| rho (vs the ic file, round trip rho,eint -> p,T -> rho) | **6.30e-09** | r = 2.40495e11 |
| p (vs the column's total p, in the untapered part rho > 1.93e-9, 3216 nodes) | **6.50e-07** | r <= 2.3396e11 |
| p (vs `p_col - (1-w) aT^4/3` over the whole domain) | 5.7e-03 | r = 2.40316e11 |

The last row is not a code error: above the taper the run's pressure is
`p_col - (1-w) aT^4/3`, a difference of two nearly equal numbers (168 vs 7563 dyn/cm^2 at
the top), so the check there is dominated by interpolating the reference column, not by
the reader.  Both gate numbers (< 1e-4) are met with four orders to spare.

Independent confirmation that the structure is the star: the run reports `tau = 2/3 at
r = 2.37122e11 cm, T = 49199 K` against the model's `R = 2.3717e11 cm, Teff = 48978 K`
(0.02 % in radius, 0.45 % in temperature), and the deep anchor at
`p = 8.35805e6 dyn/cm^2, T = 237659 K`.

## D(ii). Five turnovers in 1-D  -- FAIL

`sub_gA.sh` -> `gA.log`.  The run dies at **cycle 5, t = 56.6 s**: dt collapses to 1e-27,
the internal energy goes negative in the taper layer (`rt_use_cons gave a non-positive
internal energy in 30 cell reads`), and the top cells reach `v ~ 8e7 cm/s` (13 c_s) from
rest in 39 s.  No gate beyond (i) can be read off this run.

`sub_diag.sh` -> `diag.log`, seven 300-cycle arms that isolate it:

| arm | switch | outcome at cycle 300 |
| --- | --- | --- |
| `d_base` | -- | dead at cycle 5, t = 56.6 s, dt 6e-22 |
| `d_wall` | `outer_bc=wall` | dead at cycle 5 -- **not the outer boundary** |
| `d_nowb` | `wellbalance_dynamic=false` | dead at cycle 5 -- **not the WB scheme** |
| `d_nostrang` | `rt_strang=false` | t = 12.9 s, dt 1.9e-07 -- worse |
| `d_noforce` | `rt_rad_force=false` | **alive**, t = 556 s, dt 2.8e-02 (still degrading) |
| `d_norad` | `rt_grey=false` | **alive and healthy**, t = 2951 s |
| `d_nosponge` | -- | typo in the arm (`problem/sponge` is not a parameter); not run |

`sub_fv.sh` -> `fv.log`, `problem/rt_force_verbose=200` for 3 cycles, prints the cell
balance across the taper (i = 89..94, rho 1.68e-9 .. 5.98e-10, w 0.97 .. 0.05).  Read it
with care: the print normalises by `rt_force_grav = GM/r_in^2 = 29733 cm/s^2`, which is
**4.03x the local g** at these radii (7382 cm/s^2), and its `resid` subtracts that same
constant instead of the local g.  Converted to the local g at i = 91:
`a_p = 1.58 g`, `a_f = -0.575 g`, `a_p + a_f - g = -0.013 g`.  **The initial state is
hydrostatic to ~1 %**, so the IC is not the problem.

The remaining suspect is the mode-3 column solve on this spherical column -- see D(iii).

## D(iii). Mode 3 vs `rt_implicit_column = 0`  -- FAIL, and it localises the bug

`sub_g34.sh` -> `g34.log`, 500 cycles each, `rt_rad_force=false` (the only setting in
which mode 3 survives 500 cycles at all), `rt_col3_skip_sweep` at its default `false`:

| solver | t at cycle 500 | dt | emergent `L_rad,out / L` |
| --- | --- | --- | --- |
| `rt_implicit_column = 3` (pcr, warm 1, mixed 2) | 562.7 s | 2.1e-02 | **8.64e+05** |
| `rt_implicit_column = 0` | 6404.9 s | 1.27e+01 | **0.950** |

Gate wanted 1 %.  Mode 0 is healthy -- 5 % from `lstar`, which is the G.5 gate -- and mode
3 is six orders out and has driven its own timestep to nothing.  Mode 3's emergent flux is
right at cycle 0 (`L_rad,out/L = 0.938`) and diverges within a handful of cycles.

This is the first time the mode-3 column solve has run on a spherical mesh: the
`Av`/`Ac`/`Vc` accessors of `RTCol3` (`src/utils/two_stream_column_implicit.hpp:531-548`)
are fed the real `pcoord->area.x1f`, `pcoord->volume` and `pcoord->dx1` by
`two_stream_rt.hpp:3011-3014`, so red_giant needs no plumbing -- but the arithmetic those
accessors enable has never been exercised.  **The next piece of work is to gate the mode-3
spherical path itself** (constant-kappa radiative-equilibrium column, `L` constant in
radius, mode 3 vs mode 0), not to change the He setup.

## D(iv). Restart  -- essentially PASS (1 ulp)

`sub_g34.sh`, mode 3, `rt_rad_force=false`: 250 cycles straight (`g4_cont`) against 200
cycles + restart + 50 (`g4_rst`).  The first 201 history lines are **bitwise identical**.
The restart writes one duplicate history line at the restart time; after dropping it, the
first resumed line agrees in 9 of 10 columns bit for bit and differs in one (the x1
momentum sum) by **3.1e-16, one ulp** -- a reduction-order difference, not a state
difference.  Because this configuration is in the degenerate mode-3 regime of D(iii), that
one ulp amplifies over the following 50 cycles.  Re-run the check once mode 3 is fixed.

## Not run

* G.5's 20-turnover mass-drift / `it_mean` / `eos_fail` gates: there is no run long enough
  to read them off.
