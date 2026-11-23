# The clean M1 radiation scheme in AthenaK

This directory documents the grey radiation-hydrodynamics scheme of `src/rad_m1/` ("clean
M1", staggered layout): two-moment radiation transport (E on cells, F on faces) closed by a
variable Eddington tensor from short-characteristic ray sweeps, solved implicitly and
coupled to AthenaK's hydrodynamics. It is the scheme of the 2026 massive-star envelope
runs.

| document | for |
|---|---|
| [SCHEME.md](SCHEME.md) | what the scheme solves and how: equations, closure, discretisation, time integration, implicit solver, conservation properties, code map, how to cite |
| [BUILD.md](BUILD.md) | building on CPU, NVIDIA (CUDA) and AMD (HIP) machines; known-good toolchains |
| [RUN.md](RUN.md) | a first run, the production star decks, restarts (read the restart rule), dt retries |
| [KEYS.md](KEYS.md) | every input key the production decks use: meaning, default, production value; experimental keys |
| [OUTPUTS.md](OUTPUTS.md) | history columns, checking energy conservation, `CLEAN WARNING` messages, the end-of-run summary |

Quick start:

```bash
git submodule update --init kokkos
tests/m1_clean/run_tests.sh --quick        # builds a CPU binary, runs the quick tests (~5 min on 16 cores)
```

Inputs: `tests/m1_clean/decks/` (small test problems) and `inputs/rad_m1_clean/`
(production star templates). Tests: `tests/m1_clean/README.md`.
