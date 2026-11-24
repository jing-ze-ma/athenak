# The clean M1 radiation scheme in AthenaK

The grey radiation-hydrodynamics scheme of `src/rad_m1/` (`<rad_m1>/scheme = clean`):
two-moment radiation transport with the energy density on cells and the flux on cell faces,
closed by a variable Eddington tensor from short-characteristic ray sweeps, integrated with
an IMEX SDIRK2 scheme (implicit radiation and gas-radiation exchange, explicit AthenaK
hydrodynamics), in a well-balanced formulation for stellar envelopes.

Properties (wording, mechanism and the test of each: [SCHEME.md](SCHEME.md)):

- energy and mass conserved to round-off (floors paid from the radiation);
- angular momentum L_z: exact on the cubed sphere; on spherical polar exact except where the
  face-flux cone cut acts (booked; a fix is in progress);
- linear momentum: not conserved in general (cone cut, well-balancing residual, floor velocity
  rules; not booked), and not a discrete invariant on spherical-polar / cubed-sphere meshes;
- density positive by construction; radiation energy by the implicit operator and limiters
  (counted backstops); gas internal energy not by construction (the step is redone, or the
  radiation pays);
- well-balanced: equilibrium envelopes stay at rest to about 1e-9 of the sound speed;
- second order in space and time (asymptotically);
- bitwise restarts and dt retries;
- bitwise across MPI rank counts, and across tile layouts in tile mode (one MeshBlock along
  x1 and x2);
- CPU and GPU agree to round-off (not bitwise);
- no silent fallbacks: unsupported options are FATAL, handled events print `CLEAN WARNING`;
- Cartesian, spherical-polar and cubed-sphere meshes with one implementation.

| document | for |
|---|---|
| [SCHEME.md](SCHEME.md) | the method and its properties, code map, how to cite |
| [BUILD.md](BUILD.md) | building on CPU, NVIDIA (CUDA) and AMD (HIP); known-good toolchains |
| [RUN.md](RUN.md) | a first run, the production star decks, restarts, dt retries |
| [KEYS.md](KEYS.md) | the input keys a user sets; the restart rule |
| [OUTPUTS.md](OUTPUTS.md) | history columns, checking energy conservation, `CLEAN WARNING` messages, the end-of-run summary |
| [tests/m1_clean](../../tests/m1_clean/README.md) | the test suite |

Quick start:

```bash
git submodule update --init kokkos
tests/m1_clean/run_tests.sh --quick        # builds a CPU binary, runs the quick tests (~5 min on 16 cores)
```

Inputs: `tests/m1_clean/decks/` (small test problems) and `inputs/rad_m1_clean/`
(production star templates).
