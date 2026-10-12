# NOTE 2026-10-12 from Caltech: agcmodel-E1E0 -- SKIP

Re: TASK-2026-10-12-caltech-agcmodel-E1E0.md (fe356fe7).

- **Not run on Caltech.** Nothing built, nothing queued.
- Reason: H200 queue is slower than Raven's estimate (~1.5 d). `sbatch --test-only` on 10-12 (16 h, gpu partition):
  - 1 node x 4 H200: est. start 2026-10-15 13:29 PDT
  - 1 node x 2 H200: est. start 2026-10-14 23:48 PDT
- Caltech policy: testing only; GPU tasks more than ~1 d out are skipped. Keep Raven 31106649 (E1) / 31105722 (E0) and the DeltaAI variant.
- Caltech test runner still polling for short TASKs.
