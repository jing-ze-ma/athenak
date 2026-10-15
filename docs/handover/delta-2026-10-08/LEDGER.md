# Delta GPU-hour ledger (TASK-2026-10-08-delta-test-runner)

Account bivj-delta-gpu, 100 GPU-h deposited. Budgets: spb2 <= 40, beam <= 40, reserve >= 15.
GPU-h = elapsed x 4 GPUs (exclusive 1-node gpuA100x4 jobs; billing 1 SU per A100-hour).

| job | branch | batch / arm | submit | start | elapsed | GPU-h | note |
|---|---|---|---|---|---|---|---|
| 22759294 | rt-integration | bring-up smoke A+B | 10-08 15:00:55 | 15:01:14 | 00:01:29 | 0.10 | |
| 22760430 | rt-integration | smoke B, host -ffp-contract=off | 10-08 15:56:57 | 16:23:44 | 00:00:53 | 0.06 | bitwise = Raven 31004354 |

Totals: rt-integration 0.16, spb2 0, beam 0.
