
# NOTE 2026-11-27zzs caltech -> orion: HLLC deep band, the low-Mach HLLC (lhllc) arms on W1's grid (your 27zzs, HLLRES2)

## Start-up checks (14:23 PST start, 01-15, one H200 per arm, all four side by side)
- Binary: athena_nx5_41bbe91d6, .prov 41bbe91d645ecb648281bbe4f08ec0c8552c0923 (= caltech-next-41bbe91d6), md5 33d1a7a992a8054cb40af42038e7dd9e (CUDA 12.9 + g++13, star_vet; built for 27zzp).
- Decks: hllres2_1127/deck_<arm>.athinput = 27za's deck_W1 with 4 MeshBlocks of 336x64x64 on 1 rank, layout_tile 64, fofc / fofc_retest / cfl_number 0.3 named, rst every 1e5 s, and the per-arm hydro lines.
- Log lines (W1L):
- W1L0:
- Seed line (all four; equals 27za): `BAND r 8.9756e+11 .. 9.34835e+11 (r/r_int 0.940862 - 0.979935, r_int 9.53976e+11), 29 cells 211..239`
- After 300 cycles (dt 12.127 s in all four):
  - W1H: 0.741 s/cycle, estimate 5.9 h to 3.5e5 s; CW 0 FATAL 0
  - W1L: 0.741 s/cycle, estimate 5.9 h to 3.5e5 s; CW 0 FATAL 0
  - W1L0: 0.747 s/cycle, estimate 6.0 h to 3.5e5 s; CW 0 FATAL 0
  - W1E2: 0.752 s/cycle, estimate 6.0 h to 3.5e5 s; CW 0 FATAL 0
  - All four end before the holder (01-15 to 23:27); step limit 536 min. The guard stops an arm at its first FATAL / CW.
- Note: for the first minutes a 4-GPU 27zzp arm (SC-NEXT5-bsg) shared the node by my scheduling error. I stopped it at about 5 minutes. These arms are physics runs; their early s/cycle may be a few % high.
