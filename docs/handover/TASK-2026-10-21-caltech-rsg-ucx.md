# TASK 2026-10-21 raven -> caltech: RSG E nondeterminism + IB fault — test UCX settings (USER: yes)

Your rdiff 5140304 (RSG E deck on 4de8168d, 1 x 4 H200) shows R != C AND R2 != R (run-to-run nondeterminism), while the same checks are
BITWISE on Raven A100. Together with the link-1 IB fault (UCX transport retry exceeded, RDMA of a 5.8 MB GPU buffer in the vet_gd cs
halo), this looks like the UCX GPU-direct / rendezvous problem we found on Viper MI300A on 10-19 (large device messages corrupted;
fix there: UCX_RNDV_THRESH=inf).
1. Rerun the rdiff A case (deck E N=32, C, R, R2) twice more, 1 x 4 H200, debug QoS if possible:
   (a) export UCX_RNDV_THRESH=inf          (b) export UCX_IB_GPU_DIRECT_RDMA=no
   (c) control: unchanged env, a second time (to see if C itself is reproducible: C vs C2).
2. Report bitwise yes/no for C==C2, R==C, R2==R per variant, and the cost per cycle of each variant.
3. If (a) or (b) gives bitwise, put that export in the RSG E link script BEFORE link 1 (5126655) starts and NOTE it (hold the link
   if it would start before the test finishes). Same env for the accretor chain.
