# TASK 2026-10-22 caltech -> viper (USER): GPU memory -- breakdown, then host-side dt_retry / posfix checkpoints

Why: RSG E / HG global no longer fit 80 GB H100s (the only free GPUs on Caltech today), so they wait days for H200 nodes.
Measured (nvidia-smi / logs):
  RSG E 8 x H200, 6.75 M cells/GPU, old binary 4de8168d (10-21):            52.8 GB/GPU
  HG global 12 x H200, 5.75 M cells/GPU, old binary (10-21):                  66 GB/GPU
  HG global 4 x H200, 17 M cells/GPU, ff1f52bc + harmonic+extrap2+accept_fallback+posfix=4 (now): 131 GB/GPU
  RSG E 8 x H100 (80 GB), 6.75 M cells/GPU, ff1f52bc + same keys: OOM at startup ("failed to allocate 310.7 MiB", label dtr:hyd.w0)
  HG 6 x H100, 11.5 M cells/GPU: OOM (label m1csdxf.x2f)
Already on: vet_gd_band_ring (compact 13.1 GB + ragged 13.7 GB on rank 0 at 4 ranks RSG; dense would be 122.6 GB), vet_gd_twin_lowmem,
vet_gd_halo_compact = 2.
USER asks:
1. Breakdown: per-label device allocation of RSG E (deck rsg_E_1020_clip_harm + the 4 keys above) at 3 and 6 MeshBlocks per GPU
   (KokkosP memory-usage / kp_memory_events or your own tally), old (4de8168d) vs ff1f52bc + keys -> which arrays the fix added.
2. Cut (opt-in keys, default bitwise): keep the dt_retry checkpoint (dtr:*) in pinned HOST memory (only read on a rejected step), and
   keep/share the time2_posfix start copies (t2gsv, t2csv) with it or on the host; anything else the breakdown shows is a full-state copy
   used rarely. Gates: bitwise key-off, C2==C, R==C with the keys on, s/cycle cost on a retry-heavy stretch (RSG E 125-130 d, E2).
Goal: RSG E on 8 x H100 80 GB (3 MeshBlocks/GPU) and HG global on 12 x H100.
