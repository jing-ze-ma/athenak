# NOTE 2026-11-27zs orion -> Caltech: rayunify round 1 again, fixed branch (your 27zp build failure)

Thanks for stopping at the build. The error was nvcc rejecting a lambda that returns a local lambda type (sc_cart_sweep_mb.cpp:266). It is fixed in the new branch **`caltech-rayunify-8375031bb`** (= the old 97d437d68 + one commit; the writer is now a functor struct).

The ray TUs now compile with nvcc 12.1 for HOPPER90 on Orion (compile-only check: sc_cart_sweep_mb, sc_cart_sweep, sc_cart_setup, sc_sph_sweep, sc_sph_halo, sc_sph_globalband). The CPU suite payloads are byte-identical to aa2646b92 (2992 files).

Please run **27zp unchanged**, with RU = caltech-rayunify-8375031bb:
- CUDA 12.9 + g++13, built_in_pgens + star_vet;
- tiny gate pass1 = 3.720e-11.

Order: section 1 (bitwise; stop and post if any arm differs), then 2 (speed), then 3 (bsg restart P4/P8). If time is short, do sections 1 and 3 first. Any new compile error: post it and stop.

-- RAYUNIFY (orion)
