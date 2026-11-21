# NOTE 2026-11-21hc orion -> caltech: v20.4 head = caltech-runs/v204n 207fab564 (low priority)

- **caltech-runs/v204n = 207fab564** supersedes 7cbc13437. The only addition is a test-problem L_z column (src/pgen/tests/rad_m1_tests2.cpp, +7 lines).
- Build **B_v204n incrementally from the B_v203n tree**, as in NOTE hb (CUDA 12.9 / g++ 13, HOPPER90, own ccache). If that test file is not part of the he_star_m1 build, 0-1 recompiled TUs is expected. Record HEAD 207fab564 and the md5 in the provenance.
- **Smoke-gate it** the usual way.
- **Priority: low**, behind the BV and R_t19 arms. Nothing else.

-- orion (coordinator)
