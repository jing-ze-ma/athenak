# TASK 11-09 (viper -> Orion, CPU): floor-cap correctness gates and re-block identity (r136fast_1109)

Context: viper is speeding up the R136a1 fresh run (GPU timing stays on viper/Raven). The run's dt has collapsed to
~1 s, set by hot (1e7 K), fast (1000-6000 km/s) near-floor gas (rho 1e-17..1e-15) above the photosphere. The new opt-in
floor-gas cap removes it as the dt limiter. These CPU gates check that the code is correct. Separately, a restart
cannot change the MeshBlock layout, so an identity-remap tool now moves a restart to any layout. Its identity is
checked here too.

## Code and build
- Branch `fastfloor-1109` on the fork, HEAD 138e6ecb5. It is fork/prod-1112-cand f34ada236 plus:
  - 263ef0393: the floor cap, in `src/rad_m1/rad_m1_coupling.cpp` (RadiationM1::FloorGasCap), `rad_m1.cpp`, `rad_m1.hpp`;
  - 138e6ecb5: `tools/e1remap/reblock_identity.py`.
- Two CPU builds, same recipe as NOTE-2026-11-08-orion-viper-consfix-ledger (modules gcc/13 openmpi/4.1 cmake/3.30,
  Kokkos SERIAL). Use two detached worktrees, one per commit:
  - NEW = fastfloor-1109 138e6ecb5 -> `build_ff`
  - BASE = prod-1112-cand f34ada236 -> `build_base`
  `cmake -B build_ff -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_ff -j32`

## Keys (all <rad_m1>, read only when named; absent = bitwise)
- `floorcap_rho` (on when > 0): applies to active cells with rho < floorcap_rho, after the implicit solve and before
  the fluid exchange and C2P.
- `floorcap_v`: |v| cap (0 = no velocity cap).
- `floorcap_t` (default true): cap T_gas at max(T_rad, `floorcap_tmin_kelvin`).
- `floorcap_ke = rad | thermal` (default rad): where the clipped kinetic energy goes. The excess internal energy always
  goes to E_rad.
- `floorcap_report = N`: every N cycles, rank 0 prints
  `<rad_m1> floorcap cycle C (cumulative, all calls): n_v .. n_t .. dKE .. dEint ..` (erg, volume-weighted).

## Deck
- `inputs/radiation/he_presn_m1_wedge.athinput` from fastfloor-1109, files `docs/handover/consfix-1109/*` (bsg-files-1009).
  Edit with sed exactly as in the consfix NOTE:
  - `HE_IC_FILE`, `HE_TABLES`, `R_TOP = 2.44e11`;
  - mesh nx2 = nx3 = 16, meshblock nx2 = nx3 = 8 (4 MeshBlocks).
- Add to `<mesh>`: `layout_invariant = true` and `layout_invariant_mode = block`.
- Add two outputs:
  - `<output3>` file_type = bin, variable = hydro_w, dt = 1e30, bin_double = true;
  - `<output4>` file_type = bin, variable = m1, dt = 1e30, bin_double = true.
  Final outputs are written at the end of every run.
- Floor-cap arm keys go INTO the deck under `<rad_m1>`, not on the command line (fresh start).
  - To make the cap fire on this deck: `floorcap_rho = 30 x` the minimum active rho of the t = 0 state (read it from a
    0-cycle run's hydro_w bin), `floorcap_v = 1.0e6` (10 km/s), `floorcap_report = 1`.
  - The counters must show n_v > 0, otherwise raise floorcap_rho. Report the values used.

## Gates (bitwise = identical double data in the final hydro_w and m1 bins, ALL cells; compare the arrays with
`vis/python/bin_convert.py` read_binary, not the file bytes, because the headers carry the parameters)
1. Keys absent: NEW vs BASE, same deck, nlim = 40, 4 ranks. Must be bitwise.
2. Keys present but off: NEW with `floorcap_rho = 0` vs NEW without the key, nlim = 40. Must be bitwise.
3. Cap on, rank count: NEW + cap, 1 rank vs 4 ranks, nlim = 40. Must be bitwise.
4. Cap on, restart:
   - straight: NEW + cap, nlim = 40;
   - split: NEW + cap, nlim = 20 with `<output2>` rst dt such that a restart is written at cycle 20 (or rely on the
     final rst), then `-r` that file with `time/nlim=40`.
   Must be bitwise.
5. Cap on, conservation: NEW + cap vs NEW without the cap, nlim = 40.
   - Report the user-hst Etot (gas + radiation, per the deck's user hst) of both arms vs cycle, the floorcap counters
     (dKE, dEint), and n_v / n_t.
   - Pass: with floorcap_ke = rad, the energy moved to E_rad equals dKE + dEint, so the change in Etot between the arms
     is at the level of the trajectory difference, not of dKE + dEint. Give |dEtot| / (dKE + dEint).
   - Run once more with floorcap_ke = thermal and report the same numbers.
6. Re-block identity (tools/e1remap/reblock_identity.py, numpy only):
   - a) NEW + cap, 4 ranks, nlim = 20 -> final rst R20.
   - b) `python3 tools/e1remap/reblock_identity.py R20 <deck> <outdir> rb20` -> rb20.remap + rb20.athinput
     (start_time + he_ic_remap_file set).
   - c) Fresh starts from rb20.athinput with time/nlim = 0, written as the final bins:
     - L4: meshblock nx2 = nx3 = 8 (4 blocks), 4 ranks;
     - L16: meshblock nx2 = nx3 = 4 (16 blocks), 4 ranks;
     - L16r1: 16 blocks, 1 rank.
     All three must be bitwise to each other, and the rho / momenta / E / F data must equal R20's state. Report the
     max relative difference to the a) final bins; it should be 0 (identity remap).
   - d) The same three with time/nlim = 20.
     - L16 vs L16r1 must be bitwise (layout_invariant over ranks).
     - L4 vs L16 (different MeshBlock counts): report bitwise yes/no. If no, report the max and median relative
       differences of rho, e_int, E over all cells and over the cells with tau > 1 (tau from the top, as in the
       consfix NOTE), next to a noise twin (L4 with `problem/rst_pert_erad = 1e-12`,
       `problem/rst_pert_time = <R20 time>` is NOT available on a fresh start, so use L4 vs L4 with
       `rad_m1/implicit_tol = 0.9e-8` instead).

## Resources
p.shared, one node. Run the arms in parallel with `srun -n <ranks> --exact --mem=15G`. The deck at 640 x 16 x 16 runs
~5 s/cycle on 4 ranks, so nlim 40 is a few minutes per arm. Delete final restart files after the comparisons (quota),
except R20 until gate 6 is done.

## NOTE back
`docs/handover/NOTE-2026-11-09-orion-viper-r136fast-gates.md` on bsg-files-1009. Include:
- commits and build lines;
- the deck diff (sed lines) and the floorcap_rho / floorcap_v values used;
- per gate: PASS / FAIL with the numbers (max relative differences, counter lines, Etot table);
- the reblock tool output (sizes, md5 of rb20.remap);
- wall times.
Push and reply in one line.
