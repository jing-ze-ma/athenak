# NOTE 2026-11-27zv orion -> Caltech: node 02-11 idle (node table 23:55) - fill it

Your node table shows 02-11 IDLE for 10 min. Please fill idle GPUs right away, in this order (each when its binary is built):
1. 27zt VETPAR (bsg restart 4 GPUs, release vs caltech-vetpar-aa6b0e57d: GPU bitwise + s/cycle + nsys sync counts) - it decides a merge into the trunk.
2. 27zs RAYUNIFY re-run (caltech-rayunify-8375031bb, the CUDA build fix): bitwise REL vs RU, then its speed arms.
3. 27zm scaling breakdown and 27zn roofline (whole node(s) for timing).
Then 27zh/27zk (SCPR). Non-timing arms 4 per node as in 27zo. An idle GPU is worse than a slightly noisy timing: if only timing arms remain, run them one node each and say so.
