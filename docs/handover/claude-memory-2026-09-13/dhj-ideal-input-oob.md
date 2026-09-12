---
name: dhj-ideal-input-oob
description: "FIXED 0d1f6f6a: deep_hot_jupiter_rt_ideal_xe.athinput segfaulted in the RT within 20 cycles because wtemp is only allocated under a general EOS. Plus a one-past-the-end write affecting both inputs."
metadata:
  node_type: memory
  type: project
---

Found 2026-08-22 while timing general vs ideal EOS. **The shipped
`inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput` segfaulted in a Release build after a
handful of cycles**, inside `picket_fence_two_stream_RT`. Fixed in **0d1f6f6a**.

1. **`Hydro/MHD::wtemp` is allocated ONLY under a general EOS** (`mhd.cpp` guards the
   realloc on `IsGeneral()`), because an ideal gas has nothing to cache. The pgen read
   `wtemp_(m,k,j,i)` unconditionally at SEVEN sites as the warm start for
   `PresTempFromEint` -- indexing an EMPTY View. `PresTempFromEint` discards the guess on
   the ideal branch, so the value never mattered; the READ was the bug. Now routed through
   a `TGuess()` helper returning -1 when `wt.extent(0) == 0`.
2. **`get_wb_eos_arr` / `get_init_eos_arr`**: `for (n=0; n<N; n++) logparr(n+1) = ...`
   writes one element past an N-element array. Affects BOTH inputs, always has.

**The lesson about how to find this class of bug here.** The Release build segfaulted;
the RelWithDebInfo build of the SAME code did not -- a one-past-the-end write is silent
until it happens to land somewhere unmapped, so its visibility depends on codegen. Do not
conclude "fixed" from a build that stops crashing. The tool that found both, instantly,
was **`-D CMAKE_BUILD_TYPE=Debug`, which turns on Kokkos View bounds checking** and names
the View and the index. AddressSanitizer would not configure cleanly in this tree; Debug
did, and it is what to reach for.

Neither fix changes any answer -- the general path is bitwise identical over 100 cycles,
and both inputs are now clean under bounds checking. See [[general-eos-optimization]].
