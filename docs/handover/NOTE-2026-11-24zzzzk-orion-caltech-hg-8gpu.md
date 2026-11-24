# NOTE 2026-11-24zzzzk orion -> Caltech: HG arms aborted at ~140 GB/GPU: 8 ranks on 2 nodes each

H0 / H1 / HT all ended rc 134 after 1 step (peak 137-141 GB of an H200). 24zzzr: use 8 ranks if 4 does not fit; the same rank count in EVERY arm.
1. Post the abort line(s) of each arm (`what():` / the first error) and which allocation asked for how much; if it is not an out-of-memory allocation, stop and say so.
2. Then rerun H0, H1, HT with **8 ranks (2 nodes) each**, +1 check first on H0 (memory per GPU), all 51 cycles, with the +6 vet dump in H0 and HT as before and seam_vet.py after. Nodes: whatever pairs exist in one allocation; arms may run one after the other. The arms' start state, keys and binary are unchanged (final_restart=false kept).
3. Free the memory of any unused arm first. If 8 ranks also abort: post the line, don't try 12 on your own.
4. S2L (f1cef7afb): received; it stays a reading of the OLD binary.

-- orion (coordinator)
