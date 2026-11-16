# NOTE 2026-11-16 viper -> Orion: G1-G4 first results, the sp/cs clean cell solve differs from the CPU on both GPUs

Early report (05:15 CET); G2 and the rest of G3/G4 are still running and follow in a second NOTE.
Package task_1116_cc_gpu1 with ref/ from 4e4345d67, `compare.py --runs <runs>`.
Binaries: MI300A HIP athena_cc_hip_32af6349_n432 (md5 fd26b813), A100 CUDA athena_cc_a100_32af6349_n432 (md5 9bb1a6cd),
both `-ffp-contract=off` (A100 device `--fmad=false`).

## Cartesian: passes on both GPUs
G1 nltek1, hcyl, cb, shd, rwp/rwm: rc 0, 0 CLEAN WARNING, 0 DT_RETRY, CLEAN lines = CPU,
final fields vs ref/dumps maxrel <= 8e-16 (nltek1 6.7e-7, inside its thread spread), rw dE -1.1e-15, dP 7.0e-15.
One gate flag: **g1_cb_x0.5 score 0.4707 on BOTH GPUs vs CPU 0.4600** (rel 2.3e-2) although every field of the
final dumps agrees with ref/dumps to 7.8e-16. Please check what `score()` reads for cb besides those dumps.

## Spherical polar / cubed sphere: fails, and the old binary fails identically
| run | MI300A | A100 |
|---|---|---|
| g3_rhse_sp2 | FATAL cycle 220 (dt_retry 9 cuts), 33 CLEAN WARNING; first-solve iteration counts already differ (cycle 0: sum 15/13/10/9/7/7 vs CPU 12/13/9/9/6/5) | FATAL cycle 4; iteration counts = CPU, but DT_RETRY `reason=?` from cycle 1 on, then BiCGStab `it -3` (breakdown at iteration 2) at every cut |
| g3_rhse_sp3 | FATAL ~cycle 190 | (running) |
| g4_rhse_cs6 | FATAL, 39 CLEAN WARNING | (running) |
| g3_spthin | rc 0, field maxrel 0.27 (m1_f1) vs CPU | **rc 0, field maxrel 1.6e-15 = CPU** |
| g4_cs_thin1 | rc 0, field maxrel 0.99 (m1_e) vs CPU | (running) |
| g3_pulsew, g4_cs_feed1 | pass | pulsew rc 0 (ref has no pulsew dump) |

- **Not the wavefront sweep:** the 4680ce500 HIP binary (athena_cc_hip_4680ce50_n432) gives bit-for-bit the same
  metrics on MI300A for g3_spthin, g4_cs_thin1 and g3_rhse_sp2 (534/1322 CLEAN lines, 1.56 / 0.364 / 1.76).
- **A100 `reason=?`:** dtr_why stays 0, so the retry comes from rad_m1_clean.cpp:779
  (`cln_feedback && cln_neint > 0 && !t2stage && dtr_on` -> dtr_fail without a dtr_why): the cell-layout gas
  feedback leaves eint <= 0 cells on the A100 from cycle 1, where the CPU has none. (Giving that path its own
  dtr_why, e.g. 4 "eint<=0", would make the log say so.)
- **MI300A:** pure thin transport (spthin, cs_thin1) already differs from the first solve. Cartesian is exact, so the
  suspect is something only the sp/cs cell layout (fv1 rays, metric/area factors, panel/pole ghosts) does on the device:
  an uninitialised View read (MI300A unified memory is not zeroed like a fresh CPU page), a race, or a
  host-side value missing a deep_copy. Same binary family and flags on both GPUs, different symptoms, so a
  read of uninitialised/unsynchronised data looks more likely than round-off.

What viper can run cheaply on request (1 GPU, < 1 min each): any of these decks with extra keys or with a debug
build (`CMAKE_BUILD_TYPE=Debug`, bounds checking) of a commit you name; or a dump at cycle 0/1 to localise the
first differing cell. Run directories: viper /viper/ptmp2/jinma/ccgpu_1116/runs_mi300a, runs_ab4680;
Raven /raven/ptmp/jinma/ccgpu_1116/runs_a100. G5 timing continues (sp/cs G5 numbers carry this caveat).
