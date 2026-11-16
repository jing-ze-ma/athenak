# NOTE 2026-11-17 viper -> Orion (USER): bugs and races found on viper on 11-16/17 (what you have not had yet)

USER asked to tell you. These are for information; requests are marked.

## 1. Second MI300A race, in PRODUCTION hydro (not clean)
- Symptom: g6_sp1x1wb_base gives two very different results run to run on MI300A, with the same binary. vcs_max lands near 1e-7 or
  near 3e-3. A100 runs are bitwise. Orion's init patch (4cacdae21) is in the binary and does not remove it.
- Checksum bisect, 3 repeats on 0678599eb with debug hashes: at cycle 15 stage 2 all inputs are bitwise (u0, w0 with ghosts, u1, x1 fluxes),
  yet the x2 fluxes `hflux_x2` (src/hydro/hydro_fluxes.cpp ~l.584) differ, and everything else follows from them.
- Cause: with nx1 = 548 the x2/x3 sweeps need 3*5*ncells1*8 = 66240 B of team scratch, which is over 64 KB. TeamScratchLevel then picks
  level 1: Kokkos HIP global scratch handed from team to team with no fence. That is the same Kokkos 4.6 defect as the x1-line race
  (NOTE 588853ca9); the x1 sweep (44 KB) stays in level 0.
- Exposed on MI300A: RSG (nx1 548), HG (704), R136 (738); not M35 (~336). A100 is unaffected.
- Fix under test on viper: par_for_outer (src/athena.hpp) issues Kokkos::memory_fence before and after the body when scr_level == 1.
  6 or more repeats, gates and cost are pending. A Kokkos-side alternative is viper /viper/ptmp2/jinma/ccrace_1116/kokkos_hip_scratch_fence.patch.
- **Request:** if you build for HIP anywhere, check your level-1 TeamPolicy users. On CPU nothing changes.

## 2. push 16: clean_rot_ring default broke restarts of rotating runs (fixed 0a03129a8)
- g7_thin1_rot, 100 cycles straight vs a restart at cycle 50 or 52: with defaults (ring live), MI300A crashed with an illegal address and
  A100 was not bitwise (2.5e-14, "restart file has no implicit_predictor increment"). push 15 and ring = off are bitwise.
- Cause: the live default was decided lazily in CleanInit (first solve), after the restart reader had counted the t2 channels
  (Time2RstNchWant uses cln_ring), so the record was misread. A named key was fine, since it is read in the constructor.
- Fix 0a03129a8 (+8/-1): the default is decided in the constructor, right after the vet_gd rotation keys. GPU re-test is running.
- Your harness has no rotating deck, which is why the CPU never saw it. **Request:** add one rotating deck with a mid-run restart to reg/.

## 3. push 16 defaults, final
- ON: clean_gas_exact, clean_pc_factor, clean_rot_ring = live (rotating sets).
- Back to opt-in:
  - clean_yres: the cs6/cs24 cycle-0 stall, also on your CPU;
  - clean_fv_joint: cs_thin1 is 38-62 % slower on GPU, and +17 % on your CPU earlier.
- clean_yres = local (fork clean-cc-yresphys-viper 5dfc5c1a0) is the physical fix:
  - the cell's exact local momentum exchange, with the Newton B term consistent with it, and the correction held fixed on 0-iteration passes;
  - cs6/cs24 0/0, WB 2-12x better, momentum and energy at round-off, R136 0/0;
  - it stays opt-in for now.
- The cs6 stall mechanism: a round-off dm flip (2e-19 <-> 0) in a static thick cell, magnified by settled()'s |dv|/max|v| with
  max|v| ~ 1e-12. That test is ill-posed for static states, in all layouts.

## 4. Production rad_m1 (not clean)
- M35MMS: a rare single-cell blow-up from time2_inc_kecap. At cycle 18159 a stage increment carried KE about 1e3 times the local radiation
  energy (9.4e31 erg). One cycle later the corner cell (first active theta cell, r = 14.3 Rsun) sat at the EOS T ceiling of 1e8 K, and dt fell
  7 -> 0.35 s. A 2-rank rerun from 1.5 d does not reproduce it (trajectory-dependent). A diagnostic is running on viper.
- RSG E: the 0/1/5 cube-corner jet came back in the rerun from 6988 d, matching the old branch within 10-25 %, because the rerun kept
  cs_seam_rho_guard = 4, cs_seam_blockend = extend, cs_vertex_3cell and cs_vertex_mono (Caltech used guard 0). Analysis:
  viper /viper/ptmp2/jinma/rc18_1110/RSG_BLOB_1116.md.
  - **Possible request (user decision pending):** a hydro-only cubed-sphere corner unit test on CPU, with a blob at a corner vs a panel
    centre, and with guard 0 or vertex_3cell false.

## 5. Build gotcha
Kokkos HIP_MALLOC_ASYNC=On (the CMake default) gives binaries that segfault in MeshBlock allocation on viper MI300A. Our trees use OFF.

## Status
- push 16 candidate: fork clean-cc-push16 0a03129a8. Your CPU regression of 0678599eb is green, and the 8e-3 RT feed change is answered
  in a1886f772. The GPU gates on 0a03129a8 are finishing.
- Speed and memory branches (R136, production 1.13 s/cycle, 51 GiB):
  - clean-cc-speed2-viper: 4.63 s/cycle bitwise;
  - clean-cc-mem2-viper: 70.0 GiB bitwise, 56.8 GiB with all opt-ins;
  - a lagged FV inflow is being made second order by extrapolation (plain lag measured order 1.22).
