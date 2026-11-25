# Input decks

A deck that a test script scores lives with that script; any other deck lives here, by the problem it sets
up (not by the physics module). Decks of the radiation (VET) test suite are in `tests/vet_clean/decks/`, decks
of the regression tests in `tst/inputs/` (run by `tst/run_test_suite.py`).

| directory | what | build with `-D PROBLEM=` |
|---|---|---|
| `stars/` | VET stars: production templates `star_bsg`, `star_m35mms`, `star_t9l62` ([RUN.md](../docs/vet_clean/RUN.md)); small 3-D examples `star3d_sp` (blue supergiant wedge; its photosphere velocities are affected by the gap of the profile at r = 3.36e12-3.65e12 cm), `star3d_cs` (red supergiant, cubed sphere) with their data ([STARS.md](../docs/vet_clean/STARS.md)) | `star_vet` |
| | `red_giant_cs`: red giant envelope with the two-stream solver | `red_giant` |
| `planets/` | hot Jupiters: general EOS (`_rt_eos`), ideal gas (`_rt_ideal_xe`), reference run (`_rt_blowup`); planet keys from literature parameters: `tools/planets/` | `deep_hot_jupiter_rt` |
| `binaries/` | RY Per mass gainer: absorbing star, resolved envelope | `ry_per_accretor` |
| `cubed_sphere_tests/` | cubed-sphere hydro and MHD tests | `cs_test` (`bvals_roles_cs`: `bvals_roles_test`) |
| `spherical_polar_tests/` | spherical-polar twins of them | `sp_test` (`_mhd_lowbeta`: `sp_lowbeta`; `bvals_roles_sp`: `bvals_roles_test`) |
| `two_stream_tests/` | two-stream and correlated-k transfer: plane-parallel, spherical, thick and thin limits | `box_convection` (`box_convection_*`), `deep_hot_jupiter_rt` (`dhj_*`), `red_giant` (the others) |
| `well_balanced_tests/` | hydrostatic columns | `wb_column` (`hse_atm2d`: `hse_atm`) |
| `eos_tests/` | linear waves with the general EOS, run by `tst/scripts` | none |
| the others | AthenaK's own decks, by module | see each deck |

Run one deck (from the repository root; the head of each deck says what to look at):

```bash
cmake -S . -B build_cs -D CMAKE_BUILD_TYPE=Release -D PROBLEM=cs_test && make -C build_cs -j 8
build_cs/src/athena -i inputs/cubed_sphere_tests/cubed_sphere_rigidrot.athinput -d run_rigidrot
```

Any key can be changed on the command line, `block/key=value` (the deck must name the key). Table and profile
paths in some decks point at one machine: replace them.
