# NOTE 2026-11-16 Orion -> viper: clean init / last-face audit (follow-up to 8566dffd2)

Patch: `docs/handover/handover_1116_clean/clean-cc-init-audit.patch`. It is `git format-patch` of 2 commits on top of
int 0e2194db2, which already contains your 073378e14 (local Orion branch clean-cc-gpufix, not pushed):

- 17ce95dc5 `rad_m1: default initialisers for every scalar member; CleanFaceVector reads x2f(j+1) only off the theta faces`
- 71be376e3 `coordinates: the spherical-polar hydro source reads no x2 flux at j+1 on a 1-D mesh`

Apply with `git am` on clean-cc (it touches only rad_m1.hpp, hydro.hpp, rad_m1_clean_closure.cpp and coordinates.cpp).

## What the patch changes

| member / read | where set | where read on sp / cs / 1-D / 2-D | fix |
|---|---|---|---|
| all 304 scalar members of rad_m1.hpp with no initialiser (bool / int / Real / arrays, 247 declaration lines) | 221 only in the RadiationM1 constructor (several under conditions); about 140 only in ImplicitInit; others only in ImplicitPartitionInit, ImplicitHaloMPIInit, TauClosureInit, VetColInit, OnePassAuto, ... | any path that skips the setter reads heap garbage, which is exactly the vet_milne class of bug | each takes `{}` (zero / false), the value a CPU's fresh pages give, so CPU results stay bitwise |
| `Hydro::wb_option` (hydro.hpp) | hydro.cpp:482, only with `use_wellbalance_dynamic` | the hydro flux kernels (hydro_fluxes.cpp:140, 233, 251) and m1_fluid.hpp:77 always pass it on; UBSan reports "load of value 0xbebebebe, not a valid WBOption" | `= isodensity` (0) |
| `CleanFaceVector`: `x2f(m,j+1)` | -- | a theta face at the high edge of a **1-D** sp mesh (the fv_geo kernel); index past the end of xx2f | a theta face reads no `t2` (`t2 = t1`, unused for dir 1), the same as your `p2` fix for phi faces |
| `SrcTermsSphericalPolarHydro`: `uflx.x2f(m,IM3,k,j+1,i)` | -- | **1-D** sp: the hydro's uflx.x2f has ncells2 = 1 rows. Kokkos bounds check on g3_pulsew: `uflx.x2f [0,3,0,1,i]` of extents `[1,5,1,1,260]`. The read hit the zero IEN row, so it was harmless but out of range | reads row j on 1-D (`dj = multi_d ? 1 : 0`); the x2 fluxes there are zero, so the result is bitwise |

Not changed: `DtRetry` (dt_retry.cpp:136) still copies one uninitialised bool, `m1.tau_any_mpi`, under ASan heap fill. The cause is that my ASan build's rad_m1.cpp.o was compiled during the header edit (a build race), so it lacked the initialiser. In 17ce95dc5, `tau_any_mpi{}` is initialised (rad_m1.hpp:3201). I did not re-run that build. It is worth one check on your side.

The other clean structs (CleanArgs, CleanLevel, CleanHalo) have every scalar set by CleanMakeArgs / CleanLevelOf, so nothing changed there.

## CPU-tested on Orion
OpenMP build, gcc 13, `-ffp-contract=off`, NMAX 432. Decks from task_1116_cc_gpu1, 3-4 cycles each.

- **Bitwise against 0e2194db2.** All bin dumps and all CLEAN lines are identical for g3_spthin, g4_cs_thin1, g3_rhse_sp2, g3_pulsew, g4_cs_feed1 and g4_rhse_cs6.
- **Kokkos bounds check** (`Kokkos_ENABLE_DEBUG_BOUNDS_CHECK`, RelWithDebInfo). 0 out-of-bounds hits on g3_spthin, g4_cs_thin1, g3_rhse_sp2, g3_rhse_sp3, g3_pulsew, g4_cs_feed1, g4_rhse_cs6 and G1 (nltek1, hcyl, cb, shd, rwp), 2-3 cycles each. g3_pulsew failed before 71be376e3 and passes after it.
- **ASan + UBSan with the whole heap pre-filled** (`max_malloc_fill_size=2^31`, `malloc_fill_byte=0xbe` and `0x01`). On the 6 sp/cs decks above, the m1 dumps and CLEAN lines are identical to the zero-heap runs. The only UB report left is the tau_any_mpi line above. Caveat: that ASan binary's RadiationM1 constructor may predate the new initialisers, as explained above.

## Not tested
- No GPU builds and no runs. Neither HIP nor CUDA was compiled with this patch.
- No full-length runs: no G1-G5 gate and no regression harness.
- Multi-rank and multi-block decks (g3_rhse_sp3_mb4, g4_rhse_cs24).

## Also: the g1_cb_x0.5 score (0.4707 GPU vs 0.4600 CPU with identical dumps)
`ana/score.py` metric() for cb also reads `run.log`, with the regex `implicit_beam(\d): .*in-plane ang=([-+0-9.eE]+) deg`. It feeds that angle into the analytic exact solution in `ana/cb.py` exact().

At exactly ±45 deg the beam edges sit on the 32x32 sub-sample points, so E_L1 jumps with any change of the parsed angle:
- beam0 = -45.0001 gives **0.47070** exactly, the GPU value;
- -44.9999 gives 0.4395;
- beam1 = 45.0001 gives 0.4907.

So the GPU run.log's beam0 line differs from the CPU's `in-plane ang=-45 deg`. Please compare those two lines in the GPU run.log. This is a scoring artefact, not physics. A fix in the analysis would be to use the deck angle, or to sub-sample at offset points.
