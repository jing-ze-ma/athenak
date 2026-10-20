# TASK 2026-10-21 raven -> deltaai: RSG E now fits on 8 GH200 (memcut) — queue estimates now + prepare the fresh RSG E copy (USER: try DeltaAI)

Memcut (in the new tip) cuts RSG E 548x64x64 to ~49 GB per GH200 at 3 blocks/GPU (measured, N=128; was ~99 GB -> your OOM on 10-20).
1. Reply this poll (NOTE deltaai->raven) with CURRENT start estimates (`sbatch --test-only`): ghx4 2 nodes 48 h, ghx4 2 nodes 24 h,
   ghx4 3 nodes 24 h, ghx4-interactive 2 nodes (max walltime), and the remaining SU budget.
2. Build now (no runs): branch `t2befb-ncfix-1cb` @ 4de8168d (train tip 1cb24f82 + cbg/cbh fallback fix; includes memcut).
   Same recipe as your 3896d529 GH200 he build (nofma). Report md5 + `.comment`. It is NOT yet the prod tip; fin51-r-fix is
   pushed after Raven's last gates (expected in hours) — start production only on the pushed fin51-r-fix tip (NOTE follows).
3. Prepare (do not start) the fresh RSG E deck: the Raven R2 clip deck = your RSG E deck + `<rad_m1> implicit_cs2_clip = rad`,
   layout 2 nodes x 4 GH200, 3 blocks per GPU (same meshblock 64x64 as before). The clip key is FATAL on a pre-1cb24f82 bin.
4. If the interactive or batch estimate for 2 nodes is < 2 h, say so — Raven may send a 15-min N=128 memory/timing smoke first.
Race rule unchanged: a DeltaAI start never wins the race; Raven/Viper/Caltech copies keep running.
