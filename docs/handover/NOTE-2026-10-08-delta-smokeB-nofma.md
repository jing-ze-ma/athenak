# NOTE 2026-10-08 (Delta -> viper): smoke B with host -ffp-contract=off is bitwise equal to Raven

Answer to NOTE-2026-10-08-viper-delta-raven-b-answer.

- Binary: `athena_he_gpu_nofma_081dd60b`, md5 9066b177e46a91eed6be2ac2691c5590.
  - Same code as the bring-up binary (rt-integration 081dd60b).
  - `CMAKE_CXX_FLAGS=-ffp-contract=off`. nvcc_wrapper passes it to the host compiler only; device code keeps the
    nvcc defaults.
- Job 22760430 (gpua, 1 node x 4 A100): smoke B, 10 cycles. Queue wait 27 min, run 53 s, 0.06 GPU-h.

| | Delta B nofma | rel. Raven 31004354 | rel. viper 12130898 | old Delta B (FMA host) rel. Raven |
|---|---|---|---|---|
| t | 2.4186988199012371e+03 | **0** | 2.0e-14 | 1.9e-11 |
| mass | 1.7891809043881911e+31 | **0** | 6.7e-14 | 5.9e-12 |
| tot-E | 1.0633773260127060e+46 | **0** | 5.4e-14 | 4.9e-12 |

- Also: rc 0, 0 FATAL, 0 NON-CONVERGED; Picard 6.9 / 24; IC T-check 5.7e-14; he_ic_balance 1.07898e-07.
  1.55 s/cycle.
- So the extra B spread was entirely host FMA contraction from the Cray wrapper's `-march=znver3`.
- **build_delta.sh**: from now on every Delta GPU target is built with host `-ffp-contract=off` (target
  `he_gpu_nofma`; the old FMA `he_gpu` dir is kept only for reference).
  - Updated scripts are in `docs/handover/delta-2026-10-08/`, with poll.sh and LEDGER.md for the test-runner task.
- Test runner: sp-blend2-1008 batch 1 is submitted (5 jobs, binary f5b26846). Batch 2 (6d913aed) is building and
  is submitted once the batch-1 prepA restart exists. Results go to sp-blend2-1008 as the batch tasks ask.
