# The clean VET radiation scheme in AthenaK

The grey radiation-hydrodynamics scheme of `src/rad_vet/` (`<rad_vet>/scheme = clean`):
two-moment radiation transport with the energy density on cells and the flux on cell faces,
closed by a variable Eddington tensor from short-characteristic ray sweeps, integrated with
an IMEX SDIRK2 scheme (implicit radiation and gas-radiation exchange, explicit AthenaK
hydrodynamics), in a well-balanced formulation for stellar envelopes.

Properties (wording, mechanism and the test of each: [SCHEME.md](SCHEME.md)):

- energy and mass conserved to round-off (floors paid from the radiation: hydro with
  `eos = general`; with an ideal gas the energy of the floors is only counted);
- angular momentum L_z: exact on the cubed sphere; on spherical polar exact except where the
  face-flux cone cut acts (booked), and there too with `clean_stag_fcone_lz` + `clean_lz_spin`;
- linear momentum: not conserved in general (cone cut, well-balancing residual, floor velocity
  rules; not booked), and not a discrete invariant on spherical-polar / cubed-sphere meshes;
- density positive by construction; radiation energy by the implicit operator and limiters
  (counted backstops); gas internal energy not by construction (the step is redone, or the
  radiation pays);
- well-balanced: equilibrium envelopes stay at rest to below 1e-8 of the sound speed;
- order of accuracy: the finest pairs of the test ladders give 1.8 (time) and 1.9-2.3 (grid +
  time); also where the gas-radiation exchange is stiff (`time2_gtime_slope`);
- bitwise restarts and dt retries;
- bitwise across MPI rank counts (the state: dumps, restart files, `hydro.hst`; the sums of
  `user.hst` are not), and across tile layouts in tile mode (one MeshBlock along x1 and x2);
- CPU and GPU agree to round-off (not bitwise);
- no silent fallbacks: unsupported options are FATAL, handled events print `CLEAN WARNING`;
- Cartesian, spherical-polar and cubed-sphere meshes with one implementation.

| document | for |
|---|---|
| [QUICKSTART.md](QUICKSTART.md) | start here: download, build (CPU, NVIDIA, AMD), check, tests, plots, 3-D stars |
| [SCHEME.md](SCHEME.md) | the method and its properties, code map, how to cite |
| [BUILD.md](BUILD.md) | building on CPU, NVIDIA (CUDA) and AMD (HIP); known-good toolchains |
| [RUN.md](RUN.md) | a first run, the production star decks, restarts, dt retries |
| [KEYS.md](KEYS.md) | the input keys a user sets; the restart rule |
| [OUTPUTS.md](OUTPUTS.md) | history columns, checking energy conservation, `CLEAN WARNING` messages, the end-of-run summary |
| [PERFORMANCE.md](PERFORMANCE.md) | wall times of the star runs, comparison with Athena++ on CPUs, where the time goes |
| [tests/vet_clean](../../tests/vet_clean/README.md) | the test suite |

Quick start (details: [QUICKSTART.md](QUICKSTART.md)):

```bash
git submodule update --init kokkos
tests/vet_clean/run_tests.sh --quick        # builds a CPU binary, runs the quick tests
```

Inputs: `tests/vet_clean/decks/` (small test problems) and `inputs/stars/`
(production star templates).
