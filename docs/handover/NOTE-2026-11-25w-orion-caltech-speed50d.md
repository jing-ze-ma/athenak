# NOTE 2026-11-25w orion -> Caltech: H200 check of the speed tip with the message-volume commit (caltech-speed50d-1d10fc373)

Thanks for the corrections to the rails. Branch: 79c1cabef + cc_split (HIP only) + new ray kernel + k+ sums + dead recompute + closure-fill volume (cs with vet_rotate_avg_ema_wb). CPU 53/53 base and tip, ring tests bitwise (star3d_cs: ranks, restart, dt_retry, poison). Build star_vet at 1d10fc373 (tiny gate 3.686e-11, .prov).
1. **HG (02-03 + 02-09, 8 ranks, hegiant.00205, HGN's command) with the T6 environment** (UCX_TLS=^xpmem UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=get_zcopy UCX_MAX_RNDV_RAILS=12): ABBA of 79c1cabef (9bebf6a3) and 1d10fc373, 26 cycles each, s/cyc as in SP50 (mean, median without first / last print); hst printed-bitwise base vs tip; then one pair with the MPI timer (mt files) and one pair with the kernel timer (26 + 6 cycles).
2. **The same ABBA with the current environment (T0)**, one pair only, for the effect of the fix in the put_ppln setting.
3. **t9l62 and m35mms** (one node each): kernel-timer pair 79c1cabef / 1d10fc373, 6 + 26 cycles, hst bitwise.
If hst differs: stop and report at once. Priority over 25p leftovers; after it idle. Stop all at 08:00 PST.

-- orion (coordinator)
