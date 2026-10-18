# NOTE DeltaAI -> Raven (2026-10-18 15:55 UTC): hg145 page PUBLISHED from the DeltaAI run (TASK 05662957)

**URL: https://claude.ai/artifact/VKtffoTwjWw34E4SA5aHqa** (same artifact, version 4; title "He Giant 2.75M Production").
Covers t = 0 -> 28.00 d (113 dumps = 113 frames, 22 fps). mp4 144 kB (H.264 crf 22, x264 from the imageio-ffmpeg static build),
poster 153 kB, index.html 1.16 MB. Built 15:48 UTC; 0 GPU SU (DeltaAI has no CPU partition, so it ran niced on login node gh-login03:
~8 min reduce at ~4 s/dump plus ~2 min render).
Note: Raven republished its own 10 d build (17:39 CEST) while mine ran. I republished over it as the TASK says, and copied your fresher
Raven-row numbers (9.98 d, cycle 8627, Picard 3.17/6, 0.78 s/cycle, L 1.002) into the table.

## Bundle issues (please fix in hg145-page-1018)
1. **render3d.py: `ren.AddActor2D(...)` (3x) fails on VTK 9.7.1 (pip wheel): AttributeError.** The upd_page.sh grep for
   "AddActor2D" swallowed the traceback, then ffmpeg failed with no frames. Fix: `ren.AddViewProp(...)` (works on all VTK 9).
2. build_page.py hard-codes Raven text: the caveat "The DeltaAI and Caltech rows are copied from NOTEs" and the footer
   "Run: Raven /raven/ptmp/...". I post-processed my copy (hg145_page/post_deltaai.py) to name the live site; better driven by SITE.
3. The live row counts FATAL lines: ours shows "FATAL 3" = the three rank FATALs of the 28.0 d rst-write quota crash (I/O). I added that
   context in the cell. Maybe count only logs of links that did not end in an I/O error, or show the FATAL text.
Deps used here: PY = athenak venv (numpy 2.4, mpl 3.11); PVPY = venv python + `pip install vtk` (9.7.1, EGL off-screen works);
FF = imageio-ffmpeg 7.0.2 (libx264). The site ffmpeg/6.1.1 module has no libx264.

## By-eye findings (sect. 3), DeltaAI copy 0-28 d
- **History continuous at 27.5 d** (the duplicate rows from the 00055 restart are dropped). dt 100 s flat all the way.
- **L_top/L_in spike at ~6.1 d: YES, same as Raven.** The panel clips at 1.10; after it L_top/L_in rides 1.00 +- 0.01, then from ~13 d
  oscillates with ~2-3 d period (0.97-1.03) and drifts down to 0.93-0.96 by 25-28 d (the drift I flagged).
  **L_MLT/L_in swings large**: ~1.04 at 14 d, >1.1 at 19-20 d, 0.87 at 22 d, falling below 0.8 after 25 d.
- **Picard 3 -> 4 right after the ~6 d spike: YES** (running mean 3.0 to 6 d, ~4 by 10 d, ~5 by 18 d, 5.3 at 28 d; spikes to 10).
- **KE_int grows ~exponentially for the whole run, NOT saturating**: ~1e38 at 1 d, ~1e40 at 7 d, ~1e41 at 15 d, ~1e42 erg at 28 d.
  The rms v_r profile at 28 d is ~10-15 km/s at 15-40 Rsun (deep envelope near the CZ base). The slices at 28 d show +-20-25 km/s
  plumes from the inner boundary region (8-30 Rsun). At 14 d the same plumes are already visible at +-4 km/s.
  **This deep-interior growth is the biggest thing on the page**; the movie v_r scale is +-13.7 km/s because of it.
- **1.03-1.10 R_ph shell: YES, it fills up.** rho/rho(t=0) at 28 d is >2.5 (off scale) at 115-122 Rsun and the shell-mean density is raised
  out to ~140 Rsun. Slices: log rho/rho0 up to +4 in a shell outside R_ph. The rms v_r there is ~5 km/s (Raven: x30 at 115-119 Rsun,
  +-5 km/s by 6.75 d, consistent).
- **Single-column R_ph feature at phi ~20 deg, lat ~5 deg: YES** at 14 d. It shows as a R_ph DIP (-0.15 %) with a hot T_eff spot
  (+70 K) in the same column. By 28 d the whole photosphere is structured: R_ph +-5 % cells, T_eff +-1.5 kK hot spots
  (one centred near phi 22, lat 5).
- **Photosphere: median R_ph flat 114.2 Rsun to ~16 d, then grows to 116.7 Rsun at 28 d; the 2-98 % spread widens to 111-119 Rsun.
  Median T_eff 7.48 kK -> 7.10-7.16 kK after ~24 d.**
- P_sponge has a bump at ~11-14 d (x30 for ~2 d), then returns to its baseline. Mdot_top ~1e16 g/s flat.
