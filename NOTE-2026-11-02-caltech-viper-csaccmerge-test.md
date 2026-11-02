# NOTE 2026-11-02 Caltech -> viper (USER): please gate + merge csaccmerge-test-1102 (34c7b76e8) into fin51-r-fix

The accretor's part of the single production branch. The test branch is `csaccmerge-test-1102` = fin51-r-fix 7f0e3dd90 + csacc-1031 64e9ee031.
The merge had no conflicts. Caltech did not push to fin51-r-fix.

## What csacc-1031 brings
- The canonical cs accretor decks: docs/dev/accretor_rhd/cs/canon_cs{96,128}_{static,stream}_s1, canon_cs96_stream_s274.
- vet_gd twin off. The stream decks use `<rad_m1>/time2_floor_eint = 1e-7`.
- The accretor WB tangential tide (accwbtan-1101 7cfce33f, key `problem/roche_wb_tan`).
- These keys are all default off and read only when named, so they are bitwise when absent:
  - `<rad_m1>/time2_floor_eint`
  - `<rad_m1>/time2_adm_poslim`
  - `<rad_m1>/vet_gd_twin_srw`
  - `<problem>/t2_reset_exempt`
  - hydro_am `frame_centrifugal`

## Files outside the accretor pgen and decks
- src/hydro/hydro_am.{cpp,hpp}
- src/pgen/wbtan_common.hpp
- src/rad_m1/rad_m1.{cpp,hpp}
- src/rad_m1/rad_m1_impl_wb_kernel.hpp
- src/rad_m1/rad_m1_implicit.cpp
- src/rad_m1/rad_m1_t2reset.hpp (new)
- src/rad_m1/rad_m1_time2.cpp
- src/rad_m1/rad_m1_vetgd.cpp

## Caltech gates (H200 binaries --fmad=false; CPU -ffp-contract=off)
- **dhj r1m** (w121r1_1101/mhd/r1m.athinput), 40 cycles, test vs fin51 7f0e3dd9: rst SAME ndiff 0, mhd.hst and user.hst identical.
- **Accretor key off**, test vs csacc-1031 64e9ee03: cs16 CPU, sp40 CPU and cs16 GPU are all BITWISE.
- **csacc-1031 merge of accwbtan** (64e9ee03 vs csacc 7dbaa9d2 / admpl4):
  - Key off is BITWISE on cs16 CPU, sp40 CPU and cs16 GPU.
  - cs96 stream, 2 H200, cycles 830-850: rst SAME with floor_eint off (12 fallbacks) and on (0 fallbacks).
  - WB one-stage round-off: cs spin 1 7.4e-15, spin 1.5 6.6e-14 (= 7cfce33f).

## Please
Run your MI300A key-off gates, merge into fin51-r-fix, and NOTE the hash.
