# NOTE 2026-11-01 viper -> Caltech: WASP-121b hydro fresh-start deck, inputs, recipe and health numbers (reply to TASK-2026-11-01-caltech-viper-w121-hfresh-deck)

Everything is in `docs/handover/w121hfresh_1101/` (32 files, 14 MB; no file over 5 MB). Nothing runs on viper.

## 1. Deck and inputs
- `deck/hfresh.athinput` (md5 ac762ea0...) is the exact deck tested, am_vector ON.
  - `deck/hfresh_off.athinput` (md5 cd25eac1...) is the OFF arm, for reference only.
  - Changes vs the MHD fresh deck `pkg/s10fresh.athinput` are marked `# hfresh:`:
    - `<mhd>` becomes `<hydro>`; rsolver lhlld becomes lhllc (lhllc_cb_sensor at its default, 1).
    - Dropped: lhlld_x1_phi_min, ohmic_resistivity, max_eta, use_rkg_sts, valfven_max_cgs, problem/bbot_gauss.
    - Added `<hydro>/am_vector = true` and **`problem/am_transport = false`**. The second is required, otherwise the pgen's YieldToPgen stops with a FATAL (the correction would be applied twice). The pgen then prints a harmless "L_abs drifts" warning.
    - output3 variable becomes `hydro_w`.
- `pkg/` is the full input package, all paths relative to the run directory. Symlink `ck` and `ic` into the run dir as `job.sh` does.
  - `ck/CK3X/` (ck table Premixed_3x_g8_11_hiT2, CE table FastChem_ck_3x_int_hiT2, cia, ray, sw_flux W121) and `ic/ic_w121_3x_f5.txt`.
  - `pkg/MD5SUMS`: all 21 listed files verified OK just now against the viper package /viper/ptmp2/jinma/w121_s10fresh_1027/pkg.
  - The EOS table is not a file. `general_eos = table` is built at start-up from the composition (3x, rad_met 0.4771, eos_h2 true).
- Grid and keys, confirmed with `-n`:
  - C32 (24 MeshBlocks of 112x16x16), nx1 112, x1max 2.1e10;
  - rk2, cfl 0.3;
  - ck_impl_every 4 / xstep 8 / once, ck_nquad 2;
  - cs_seam_flux positive, top sponge on, bottom off, etotgrav true;
  - tlim 300 rot; outputs: hst 0.01 rot, log 0.1 rot, bin 0.5 rot, rst 1 rot.
- No `<job>/recipe`, so the new defaults apply. The start-up block prints `<hydro>/am_transport = true (default)` and `<time>/dt_retry = true`.

## 2. Binary recipe (`build.sh`)
- Commit fork/fin51-r-fix **182477553**, `-D PROBLEM=deep_hot_jupiter_rt`, `-D Athena_ENABLE_MPI=ON`, `-D CMAKE_BUILD_TYPE=Release`. Our viper md5 is 38c3fa6c, for reference only.
- **Viper-only:**
  - `Kokkos_ENABLE_HIP`, `Kokkos_ARCH_AMD_GFX942_APU`, `Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF`;
  - modules gcc/16 rocm/7.2 openmpi_gpu/5.0, hipcc 7.2.4;
  - run env `HSA_XNACK=1`, `HSA_NO_SCRATCH_RECLAIM=1`.
  - Replace these with your CUDA recipe (Kokkos_ENABLE_CUDA, your Hopper arch flag, nvcc_wrapper).
- fmad-1031 did nothing on HIP. On your CUDA build with `-ffp-contract=off` it turns device FMA off unless you set `-D Athena_DEVICE_FMA=ON`.
- Later fin51-r-fix commits up to 3e0bd80a9 are rad_m1 / halo work that the dhj hydro deck does not touch. Basing on them is fine.

## 3. Health numbers, viper MI300A, 1 node x 2 GPU, t = 0 to 2.7 rot (`RESULTS.md`, `ana_out.txt`, `ledger_out.txt`)
- **Speed and dt:** 0.0277-0.0279 s/cycle (am_vector off: 0.0273). dt settles at 13.79 s (mean 13.6-13.7 s), about 8050 cycles/rot.
- **Health:**
  - no FATAL in 6 runs;
  - eos_fail, eos_vceil, c2p_it and fofc all 0;
  - floor counts per 10 min (about 2.6 rot, ON / OFF): dfloor 4.9e8 / 5.0e8, efloor 1.10e9 / 1.13e9, tfloor 3.6e8 / 3.8e8. Efloor energy is only ~2e-8 Lsw after rot 1.
- **Memory:** GTT peak 20.5 GiB per GPU.
- **Energy ledger, ON arm (units Lsw_abs = 2.79e30 erg/s), rot 0-1 / 1-2 / 2-2.6:**
  - dE/dt -0.0328 / -0.0385 / -0.0333: a fresh IC, cooling;
  - Lir_top 1.0366 / 1.0419 / 1.0367;
  - Etot_top +0.0366 / +0.0419 / +0.0367;
  - Etot_bot +3.35e-3 in all three;
  - Lsw_refl 0.383;
  - residual +3.9e-4 (rot 0-1, start-up) / +4.8e-5 / +4.2e-5;
  - mass drift -2.5e-13.
  - The OFF arm agrees to 4 digits.
- **Angular momentum:**
  - |Lrel|/(Iw Omega) ~ 2e-6; Lrel 9.25e37 -> -1.01e37 at 2.6 rot; Tcons -1.5e25.
  - The spin reservoir is nonzero in every active cell, with peak 9e34 -> 2.4e35 from rot 1 to 3.
  - ON restarts are 8.9 MB larger than OFF (the reservoir slab).

## 4. What I would change or watch for the 300-rot production
1. **Keep the deck as is for the solver keys.** The header comment block still describes the MHD B seed. It is cosmetic only; please ignore it or trim it.
2. **Spin reservoir is missing from the AM history.** The am_hst columns (Lrel, Iw, T*) do not include it. To check total-L conservation with am_vector on, sum the reservoir from the rst files (`scripts/spincheck.py` reads the slab), or add a column.
3. **Fix the FMA choice at launch.** Keep it the same for the whole chain, so restarts and later arms stay comparable.
4. **Watch the deep velocity checkerboard.** It is the hydro lhllc legacy we had to filter out of the old rot-300 state. The fresh start plus lhllc_cb_sensor = 1 should keep it away, but its regrowth source is the panel seams.
   - Measure the grid-scale content (odd-even fraction of v_h, shells 2-30, equator) at rot ~20 / 50 / 100. Our yardstick is oe < 0.3 and no rise at seams. On the old state it was 0.66-0.9.
   - Your seam A/B decides the seam keys; we would adopt its result.
5. **Roche port.** The IC (`ic_w121_3x_f5.txt`) is a 1-D spherical RCE profile. With the tide on, the start is no longer hydrostatic unless the pgen maps it onto the equipotentials as you describe. Check that the first ~0.1 rot has no large radial transient (vr rms in the deep shells against the tide-off arm).
6. **dt may still drop as the flow develops.** We measured only 2.7 rot; dt_retry is on by default.
