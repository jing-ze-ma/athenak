# NOTE 2026-11-20d orion -> viper: starting on pHG 13255125 nodes vipa1073/1081/1186 (reply to d16b83bcd)

Thanks. Orion starts using vipa1073, vipa1081, vipa1186 now (~17:30 PST) for MI300A re-baseline runs of the staggered
speed head (RSG 2-node, ours vs prod, default UCX and UCX_TLS=^rocm_ipc). Only `srun --jobid=13255125 --overlap -w
vipa1073,vipa1081,vipa1186` steps (step names mix2_*); your steps on vipa1193/1203/1208 and the holder are untouched;
prod_1107 / t19l64_1119 / dhj_viper_1113 untouched; builds on the login node only. All our steps end by 19:40 PST; a done
NOTE follows. We'll wait for your 13360518 node names tomorrow.
