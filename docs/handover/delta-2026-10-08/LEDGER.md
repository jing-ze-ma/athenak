# Delta GPU-hour ledger (TASK-2026-10-08-delta-test-runner)

Account bivj-delta-gpu, 100 GPU-h deposited. Budgets: spb2 <= 40, beam <= 40, csfloor <= 10, reserve >= 5.
GPU-h = elapsed x 4 GPUs (exclusive 1-node gpuA100x4 jobs; billing 1 SU per A100-hour).

| job | branch | batch / arm | submit | start | elapsed | GPU-h | note |
|---|---|---|---|---|---|---|---|
| 22759294 | rt-integration | bring-up smoke A+B | 10-08 15:00:55 | 15:01:14 | 00:01:29 | 0.10 | |
| 22760430 | rt-integration | smoke B, host -ffp-contract=off | 10-08 15:56:57 | 16:23:44 | 00:00:53 | 0.06 | bitwise = Raven 31004354 |
| 22760753 | sp-blend2-1008 | batch1 prepA | 10-08 16:21 | 16:24:43 | 00:00:14 | 0.02 | FATAL output1/dcycle not found |
| 22760754-57 | sp-blend2-1008 | batch1 A1-A3, B1 | 10-08 16:21 | - | 0 | 0 | cancelled after prepA failed |
| 22761186 | sp-blend2-1008 | batch1 prepA retry (Delta-fixed inputs), regular | 10-08 16:35 | - | 0 | 0 | cancelled: interactive copy started first |
| 22761303 | sp-blend2-1008 | batch1 prepA retry, interactive (2x rate) | 10-08 16:42 | 16:47:44 | 00:00:32 | 0.07 | FATAL: 2 rst blocks (input output5 + output6) |
| 22761187-90 | sp-blend2-1008 | batch1 retry A1-A3, B1 | 10-08 16:35 | - | 0 | 0 | cancelled; superseded by batch1b |
| 22761659 | sp-blend2-1008 | batch1b prepA (b6240ae4 scripts, nofma 6d913aed), interactive (2x rate) | 10-08 16:55 | 17:11:16 | 00:12:09 | 1.62 | rc=0, rst 00000-00002 written |
| 22761660-66 | sp-blend2-1008 | batch1b A1-A3, B1; batch2b A4-A5, B2 | 10-08 16:55 | - | 0 | 0 | cancelled while PENDING per NOTE-2026-10-09-viper-spb2-cancel-delta-arms (Raven answered 1b/2b); batches closed |
| 22772932 | rad-beam-1008 | batch3 BSG hrdet gate, 4 nodes x 4 A100, interactive (2x rate) | 10-09 09:08 | - | 0 | 0 | cancelled while PENDING (est. start 10-10 19:46): DeltaAI answered the gate first (bsg-files-1009 bb70ef9e, viper ack 485a12e4) |
| 22774481 | rad-beam-1008 | mem batch1 (mem_r4, base_r4, mem_c4), 1 node, interactive (2x rate) | 10-09 10:27 | - | 0 | 0 | cancelled while PENDING (est. start 11:47): DeltaAI answered first (bsg-files-1009 2be3c76a) |

Totals: rt-integration 0.16, spb2 1.71, beam 0 (batch3 + mem batch1 built, never ran), csfloor 0.
