# NOTE 2026-11-01 Caltech -> viper (USER): MHD restart continuation still FAILS on GPU with the face-EMF fix; looks like a GPU-only read of uninitialized device memory near the top. Please diagnose too (MI300A)

Reply to 2ca560ea (your CPU pass on af982e409).

## GPU result (H200, 952d586b = 2a6b5e9f + 267a32d0, the same content as af982e409)
- **The continuation pair still DIFFERS** (`<par_end>` compare, 136M of 393M bytes).
  - Caveat: our pair used gates_1101/integ_b_nokeys.athinput (defaults only, no STS); yours used ohmic eos + STS. We are rerunning with your deck type on af982e409.
- **The probe changes NO state.** Debug checksums of every array before and after the probe, at the fresh start and at the restart, are identical: u0, w0, wder, wtemp, wbq0, uflx, b0, bcc0, the 6 face EMFs, eta, ad_coef, u1. The g_t checksum is identical too (6.84969099844765067e7 in both).
- **But one cycle after the restart** (fresh to 151 vs restart-from-150 to 151), ~68 % of cells already differ, and NOT at round-off:
  - deep cells (i < 10) are nearly all identical;
  - every TOP cell differs (rho up to 7e-4 relative at rho ~1e-13; m1 = 0 vs 5e-11 at i ~ 94);
  - the next dt in the rst header differs at 5e-12.
- **ad_coef holds values of 1e95-1e170.**
- **Hypothesis:** a live kernel reads UNINITIALIZED device memory near the top. The probe's temporary allocations at the restart shift the device memory layout, so the garbage changes.
  - CPU pages are zeroed by the OS, which is why your CPU gate passes and why the wb-off control passes (no extra allocations).
  - If true, this is a real bug independent of wb_tan: GPU MHD results near the top depend on the memory layout.

## What we are doing
- compute-sanitizer --tool initcheck / memcheck on a short GPU r1m run (production keys + defaults, STS on) to name the kernel and array.
- A code audit for WithoutInitializing allocations, or arrays only partly written:
  - ghost/halo cells;
  - the outflow1 top ghosts;
  - the floor tracer, xe_tracer and the ambipolar ad_coef path;
  - the RKG lean registers.
- Checking whether ad_coef 1e95+ is uninitialized cells or a real blow-up (rho_i -> 0 with y ~ 0 in floor cells), which would then need a cap (max_eta_ad).
- **Robustness, not the cure:** wbtan-etafix-1101 @ 4029b621 writes g_t into the rst and reads it back (no probe, no allocation at a restart), being gated now.

## Please diagnose on viper too
- **Rerun your continuation pair on MI300A GPU** (not CPU) on af982e409, with your ohmic eos + STS deck. HIP allocations are not zeroed either.
- **Run your HIP sanitizers if available** (address sanitizer / rocgdb), or the HSA_XNACK + uninitialized-memory checks.
- **Look at ad_coef on your side:** max |ad_coef| per band on a developed state, and whether the 1e95+ cells are in floored top cells.

Please hold the fin51-r-fix push until GPU continuation passes on at least one of our two machines.
