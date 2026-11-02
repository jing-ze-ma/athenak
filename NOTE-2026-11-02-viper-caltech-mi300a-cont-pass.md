# NOTE 2026-11-02 viper -> Caltech (USER): GPU continuation PASSES on MI300A (af982e409); ad_coef has no uninitialised read; please run a determinism test on H200

## Viper MI300A results on af982e409
- **Fresh r1m continuation** (bare deck = all defaults, ohmic eos + use_rkg_sts = true, 2 ranks), 20 + 20 vs 40 cycles: rst `<par_end>` payload IDENTICAL.
- **The no-STS variant (your deck type)** is running now.
- **GPU 150-cycle restart gate** (s10/b100 vs d3f878be): rst data + mhd.hst IDENTICAL.
- **Raven CPU:** s5_mhd suite 41 passed.

## ad_coef review (our code)
- **Set over all cells incl. ghosts,** from ConToPrim, from TracerSource every cycle, and from the RKG stage.
- **Allocated with Kokkos::realloc (zeroed)** and deep_copied, not WithoutInitializing.
- **The 1e95-1e170 values are by design:** rho_i -> 0 gives coef up to 1e300, in cells with tiny x_e (the cold top and deep layers). Every reader uses fmin(coef*B^2, max_eta_ad), so inf caps and 1e300*0 = 0. There is no NaN path.
- **The EMF reads only cells inside the set range.**

So we find no uninitialised read there.

## Please test determinism on H200
- **Run the SAME 300-cycle straight run twice** (same binary, same deck, same GPU type) and compare the rst payload.
- **Why:** "m1 0 vs 5e-11 at the top" looks like the valfven_noinfall_fac threshold flipping. It zeroes inward m1 where rho <= 4 B^2/v_A^2, which amplifies any last-bit difference.
- **If two straight runs already differ on H200,** the cause is GPU run-to-run nondeterminism (atomics or reduction order), not restart continuation.
- On MI300A and CPU, every pair we ran was identical.
- Your sanitizer run will confirm or rule out uninitialised memory independently.

## Push
- We push af982e409 to fin51-r-fix once the MI300A no-STS continuation also passes. That meets your condition ("passes on at least one machine").
- If your determinism or sanitizer test finds a real bug, it gets a follow-up fix.
