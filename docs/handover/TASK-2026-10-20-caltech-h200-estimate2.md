# TASK 2026-10-20 raven -> caltech: fresh H200 start estimates for the global RSG / He giant (reply quickly; nothing to build yet)

USER wants the fresh global RSG (cs N=128, design E, 24 MeshBlocks 548x64x64) and the global He giant (24 blocks 704x64x64, ~33 GB each)
to start as early as possible; we race sites (first to start wins). Current: RSG Raven ~15 h, HG Viper after its maintenance (>6.5 h).
Please run `sbatch --test-only` on the H200 gpu partition, 24 h, for 8 H200 (2 nodes) and 6 H200 if your nodes allow, plus 12 h variants,
and reply in a NOTE caltech->raven with the start estimates. If one is earlier than ~10 h, say so; we'd then send a build TASK
(train tip prodmerge-1019b f9fb3406, decks staged in docs/handover).
