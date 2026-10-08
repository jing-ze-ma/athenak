# Delta GPU-hour ledger (TASK-2026-10-08-delta-test-runner)

Account bivj-delta-gpu, 100 GPU-h deposited. Budgets: spb2 <= 40, beam <= 40, reserve >= 15.
GPU-h = elapsed x 4 GPUs (exclusive 1-node gpuA100x4 jobs; billing 1 SU per A100-hour).

| job | branch | batch / arm | submit | start | elapsed | GPU-h | note |
|---|---|---|---|---|---|---|---|
| 22759294 | rt-integration | bring-up smoke A+B | 10-08 15:00:55 | 15:01:14 | 00:01:29 | 0.10 | |
| 22760430 | rt-integration | smoke B, host -ffp-contract=off | 10-08 15:56:57 | 16:23:44 | 00:00:53 | 0.06 | bitwise = Raven 31004354 |
| 22760753 | sp-blend2-1008 | batch1 prepA | 10-08 16:21 | 16:24:43 | 00:00:14 | 0.02 | FATAL output1/dcycle not found |
| 22760754-57 | sp-blend2-1008 | batch1 A1-A3, B1 | 10-08 16:21 | - | 0 | 0 | cancelled after prepA failed |

Totals: rt-integration 0.16, spb2 0.02, beam 0.
