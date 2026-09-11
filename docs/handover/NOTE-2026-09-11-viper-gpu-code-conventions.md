# GPU (HIP / MI300A) code conventions — read before writing kernels on orion

Written on viper, 2026-09-11. Two things that compile and pass on orion's CPU/CUDA builds
but break, or silently change answers, on the viper HIP build (hipcc, ROCm 6.3.4, gfx942).
Both have now bitten twice; please keep them in mind when writing code on orion.

## 1. Kokkos DualView sync: use the space-free idiom

Do NOT write

```cpp
view.template modify<DevExeSpace>();
view.template sync<HostMemSpace>();
```

The explicit-space templates fail a Kokkos `static_assert` under hipcc (the memory space
you name is not the DualView's device/host space on HIP, and the template is rejected at
compile time). Write instead

```cpp
view.modify_device();   // device copy was written
view.sync_host();       // bring it to the host
// or the other direction:
view.modify_host();
view.sync_device();
```

This is what `radiation_tetrad.cpp` already does. Hit in cdd7d2a5 (the hydro dt
diagnostic from the orion merge would not build on HIP) and again in 27ca5b13 (the
cubed-sphere seam table, written on a CPU build). A CPU build accepts both spellings,
so this is invisible until someone builds on viper.

## 2. hipcc changes round-off when arithmetic moves between kernels

Splitting one kernel into two, or moving an expression into/out of an inlined function,
can change what hipcc FMA-contracts, even when every floating-point expression is
character-for-character unchanged. The CPU build (g++) stays bitwise identical; the HIP
build differs at the last bit and the difference grows chaotically in the dhj runs.
Seen in the orion merge (caad9247, the `floors_legacy` gate exists because of it) and in
77de5618 (the implicit radial solve split: 1 ulp in one cell after 3 cycles, deterministic,
no race).

Consequences:
- "bitwise identical on CPU" does not imply bitwise identical on GPU. If a change must
  be bitwise on the production machine, check it on that machine (old-binary-twice as the
  run-to-run control, then old-vs-new).
- A GPU last-bit difference after a pure refactor is expected; a difference that is NOT
  reproducible run-to-run is a race and must be fixed.
- Restructurings that preserve the per-cell arithmetic and only hoist loop-invariant work
  (1bd31298, 27ca5b13) stayed bitwise on HIP; those that split a kernel did not.

## 3. Cheap habits that would have avoided both

- Build with `-D CMAKE_BUILD_TYPE=Debug` once (Kokkos bounds checking) before calling a
  new View index scheme done; the seam table had a ghost-range overrun only Debug caught.
- Copy everything a lambda needs into locals before `par_for`; never touch `this`.
- Keep `Kokkos::realloc` labels unique so `rocprofv3 --kokkos-trace` can name the kernel.
