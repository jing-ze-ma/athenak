# NOTE 2026-10-22 viper -> caltech (USER): GPU memory keys merged; fin51-r-fix = 0429ff46b

fin51-r-fix fast-forwarded: 1976bac68 -> f9d6d996a (posfix for am_vector ed84cfed8 + N^2 seed / auto ramp v1 cf26bc2f7) -> 0429ff46b
(+ gpumem-1022 3617c794). Ramp rule v2 (db6e6c5c) is NOT merged yet (see the sigma_max note).
New keys, all opt-in, default bitwise:
  <time>/dt_retry_host = true            dt_retry checkpoint (dtr:*) in pinned host memory (2.34-2.54 GB/MB off the GPU)
  <time>/dt_retry_host_skip_clean        default true under dt_retry_host: vra_d (1.15 GB/MB) copied only when written
  <rad_m1>/time2_posfix_host = true      posfix start copies t2gsv/t2csv in pinned host memory (0.11-0.2 GB/MB)
  <rad_m1>/restart_pack_per_mb = true    restart writer stages the hesdirk2 channels per MeshBlock (-0.8 GB/MB peak transient)
RECOMMENDED for RSG E / HG on H100: all four (needed for HG on 6xH100; optional margin for RSG E on 8xH100).
CUDA (A100) per MeshBlock (548x64x64): keys off 16.7 steady / 17.4 peak GB/MB; keys on 13.9 / 13.95 GB/MB (MI300A agrees).
Production layout measured on MI300A (N=128, 8 ranks, 3 MB/GPU): new binary 53.9 / 56.3 GB steady / peak (old 4de8168d 52.5 / 53.3).
Predicted per GPU (+~0.6 GB CUDA context): RSG E 8xH100 ~56-57 GB peak keys off (fits; your OOM = gres-less step, as you found),
~48 GB keys on. HG 6xH100 ~86-91 off (OOM) / ~72 on (tight); HG 12xH100 ~45 off / ~37 on.
What the fix added vs 4de8168d: +0.5 GB/MB steady (Time2Init extrap2/posfix arrays +0.35, posfix t2gsv +0.11), +0.5 GB/MB restart-write
transient. Big items are old: implicit work array iw 2.0, other ImplicitInit 2.0, dtr checkpoint 2.34, vet_gd ragged 2.23, Time2Init
1.28, vra_d 1.15 GB/MB. Not cut (needs design): iw (92 comps, used in every solve).
Gates: key-off bitwise vs tip; keys on == keys off bitwise on a 60-cycle RSG E N=32 stretch with 3 DT_RETRY + 3 posfix; C2==C; R==C
(MI300A + A100 + CPU). Combination 0429ff46 vs f9d6d996a: RSG E det, accretor D_B1v6, N=32 stretch bitwise.
Cost: MI300A +3.9 % s/cycle (checkpoint save 45 ms/cycle), A100 +3.1 % (2 interleaved pairs); H100 3 MB/GPU expected ~0.1-0.2 s/cycle.
Async save is not safe (stage 1 writes the copied arrays), so the save stays synchronous.
MI300A binaries (viper): /viper/ptmp2/jinma/gpumem_1022/bin/athena_he_apu_0429ff46_r71 (886995f9), athena_apu_ry3d_0429ff46_r71 (f4f29653).
