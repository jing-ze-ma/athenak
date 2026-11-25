# NOTE 2026-11-25v orion -> Caltech: the transport result, and the follow-up

Excellent test (8cfff61db): get_zcopy + 12 rails takes HG from 3.33 to 2.775 s/cyc (-17 %, bitwise, Waitall 1.1 -> 0.58 s). Keep these as the HG two-node environment for all further HG arms: `UCX_TLS=^xpmem UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=get_zcopy UCX_MAX_RNDV_RAILS=12` (T6). No more variants. (1) one arm of the OLD HG binary (HGO's 1d1972805, its own keys) with T6 for 26 cycles, bitwise to HGO where it overlaps, so that old / new are compared with the same environment. (2) one-node stars do not need it (Waitall 25 ms): do not change their launch. (3) the production stars are untouched (on any node): this is for tests.
The kernel-timer / speed tips: as soon as I post the message-volume tip (Orion gates running), the HG arms run with T6, base and tip. Stop at 08:00 PST as ordered.

-- orion (coordinator)
