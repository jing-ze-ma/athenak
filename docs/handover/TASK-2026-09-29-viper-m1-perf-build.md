# TASK for viper: compile branch m1-perf-0928 (HIP) and run one He-box check, BEFORE it is merged into rt-integration

**User 09-29:** the DeltaAI M1 speed-up branch `m1-perf-0928` goes into rt-integration only after viper has built it
and checked the He box. The branch is on the fork. It is based on 401875f0, and it has no conflicts with
rt-integration b3a3f53f, whose changes since 401875f0 are all docs. Full results:
`docs/handover/NOTE-2026-09-29-m1-perf.md` **on that branch**.

**What the branch changes, for HIP:**
- `<rad_m1>/implicit_opac_newton` is now **DEFAULT ON** wherever `implicit_opac_update` and `implicit_gas_newton`
  are both on. That includes the He box (`force_reference = wb_arad`), in Cartesian and sp rows, hydro and MHD.
  Restarts whose file lacks the key keep it off.
  - On GH200 this is -58..-65 % ms/cycle on the He box, with Picard passes 6.8 → 2.0.
  - Time order and ringing are unchanged: from an evolved restart, dens/ener order 1.9-2.0 for both, and the two
    agree to 3-4 digits.
- Also in the branch:
  - A fused Picard-residual reduction (bitwise).
  - The CUDA event fix. It is CUDA only; the HIP `M1Evt` is unchanged.
  - Opt-in `implicit_halo_ipc`: CUDA only, behind `KOKKOS_ENABLE_CUDA`, **never compiled on HIP yet**.
  - Diagnostic key `dbg_opac_part` (default 0).
  - `rad_m1_wedge`: `wg_ic = grey` accepts a non-const opacity on a restart.

## Steps (viper, ROCm 7.2, one short apudev job)

1. `git fetch` the fork. Build `box_convection` from `m1-perf-0928` with the viper box build (ROCm 7.2 target).
   **Report any HIP compile error**; the new code is in `src/rad_m1/rad_m1_{implicit,krylov}.cpp` and `.hpp`.
2. He box at 1 GPU, 200 cycles, `bench-2026-09-29-hebox` input, twice:
   - (a) the input as shipped. It does not name the key, so the new default is on.
   - (b) the same with `implicit_opac_newton = false` added under `<rad_m1>`. The command line cannot add keys.
3. Check:
   - rc 0, no FATAL, 0 NON-CONVERGED.
   - Picard mean ~2.0 in (a), ~6.8 in (b), from the `implicit transport: solves=... Picard iterations mean=` line.
   - The hst of (a) vs (b) agrees to tolerance level (DeltaAI: <= 2e-7 relative in mass, energy, x1-momentum and
     KE; the near-zero transverse momenta may differ at the noise level).
   - ms/cycle of (a) vs (b).
4. Optional: 2 GPUs, the same pair.
5. Push a one-line status to `docs/handover/` (or relay via the user): build OK or not, the Picard means, the
   ms/cycle, and the hst agreement. The DeltaAI session merges into rt-integration on the user's go-ahead after that.
