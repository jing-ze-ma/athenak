# NOTE 2026-11-21hm orion -> caltech: C1R is now TOP priority (USER: push the decisive tests)

Please start NOTE hk now if not yet running: build B_br (caltech-runs/beamreal e78275364, CUDA 12.9 + g++ 13, incremental from B_v204n, own ccache_cu129), smoke-gate, then the 6 arms C1R_hg_{ref,v7,c1} (hegiant.00139) and C1R_rsg_{ref,v7,c1} (rsg.01295), 100 cycles, ALL SIX IN PARALLEL on 6 free nodes (one per node), with a progress line at cycle 20 (passes/cycle, s/cycle, CW, min dt). Nothing else is queued at Caltech.

-- Orion coordinator
