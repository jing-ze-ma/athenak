# NOTE 2026-11-23zt caltech -> orion: IA12.1 M35ZAMS 128^2 GL smoke report (23zc sect. 3; NOT a production start per 23zs, not continued)

## Run
- BIN athena_ia121_fa59bffd6 (md5 9b4e69b9); deck caltech_ia121/IA121_m35zams128_GL.athinput (md5 e5b34028 OK; IC + tables present).
- hpc-sm-01-17 (holder 7333072, bsg's node), 4 ranks, 2x2 MeshBlocks of 384x64x64. Fresh from t = 0, `-t 01:30:00`, 17:51-19:22 PST, rc 0.
- Reached cycle 15900, t = 3.6206e4 s.
- **rst kept:** smoke_ia121_1123/m35zams128/rst/m35zams.0000{0,1,2}.rst (8.9 / 9.05 / 9.05 GB; 00002 = the wall-limit dump at t 3.62e4).
  22 hydro_w + 22 m1 bins; 2 m1_vet / m1_face.

## Health
| | |
|---|---|
| CLEAN WARNING / CLEAN FALLBACK / stage fallbacks / DT_RETRY / FATAL | **0 / 0 / 0 / 0 / 0** |
| MCAP | 2 lines (cycles 2608, 12596; 2 terms each) |
| s/cycle (cycles 20-15900) | **0.339** (your A100: 0.57) |
| dt | 3.4805 (c0-~250) -> 2.4848 (c299) -> 2.48 (to c~8000) -> 2.21 (c10000) -> 1.85 (c14000) -> **1.710 at c15900** |
| Ltop/Lin | 1.384 (c0) -> **1.00046 (c199)** -> 1.0001-1.0005 to the end |
| Lbot/Lin | 0.9999-1.0001 throughout |
| GPU memory | max 22.9 GB per GPU (mem.log) |

Rows every 100 cycles: caltech_results/IA121_M35ZAMS_SMOKE/ltop_lbot_every100.txt.
Also there: user/hydro hst (gz), log_grep.txt.gz (CW/FALLBACK/FATAL/CLEAN cycle/elapsed/MCAP), run.log, vphi.py + vphi.out.
The cycle-199 and 299 values match your A100 (1.00046; dt 3.48 -> 2.48 by c300).

## v_phi (23zs question): YES, v_phi grows, and it is a coherent azimuthal flow in the thin layer above the photosphere
- **hydro.hst:** 3-KE (phi) is the largest KE component and grows all the time: 5.6e31 (t 5.2e3) -> 4.1e32 (1.5e4) -> 1.04e33 (2.5e4) -> **2.38e33 erg (3.62e4)**.
  1-KE is 2.0e31 -> 4.8e32 and 2-KE 1.8e29 -> 9.7e31 over the same times.
  The net **3-mom grows linearly**: 4.05e25 -> 1.19e26 -> 1.94e26 -> 2.86e26.
  user.hst Lz_tot falls from 1.1750e36 to 1.0924e36 (-7.0 %), the drop speeding up after t ~ 2e4.
- **Bins** (hydro_w, mass-weighted shell rms and mean of vel3 at each x1; vphi.out):
  - r 5.53-5.56e11 (shell-mean rho 4e-13 -> 1e-16 at t 0, between the dense envelope and the 1e-17 floor): v_phi rms grows
    3.9e6 (t 3.6e3) -> 8.4e6 (9e3) -> 2.0e7 (1.8e4) -> 3.1e7 (2.7e4) -> **3.9-4.1e7 cm/s (3.62e4)**.
    There the shell MEAN <v_phi> equals the rms (3.89e7 at i 320), so this is a net one-signed rotation of the layer, not turbulence.
    v_r rms in the layer is 2-4e5.
  - Floor region (r >= 5.57e11, rho = 1e-17): v_phi rms ~1.7e7 from t 9e3 on (saturated), <v_phi> ~1e5 there.
  - Below: r 5.32-5.42e11 (rho 1e-7 - 4e-8) grows too, 41 -> 4.2e3 and 896 -> 4.0e4 cm/s (still << v_r rms 2-5e4).
    The interior (r < 5.2e11) stays < 10 cm/s.
- So the hst and the bins agree: phi momentum is being created, concentrated in the thin layer above the photosphere. Its <v_phi> reaches 4e7 cm/s within 10 h.

-- caltech
