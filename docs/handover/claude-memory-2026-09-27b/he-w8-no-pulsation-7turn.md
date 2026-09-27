---
name: he-w8-no-pulsation-7turn
description: He box_w8 at 7.7 turnovers has NO pulsation - only the 148 s box cavity mode, excited by the initial transient and killed by the bottom sponge (e-fold 600 s); rt_surface cadence 75 s cannot resolve 150 s
metadata:
  type: project
---

analysis_0916/pulse_he_w8 (09-17 05:10): k=0 spectra of mean top flux, <v_z>=mom1/mass, KE1 in
box_w8 to 7.7 turnovers. Only coherent line: P = 148.3 s in <v_z> (coherence 980), amplitude
4e-5 c_s in turnover 0-1 then 1e-12..1e-10 c_s (noise floor) from turnover 2 on, e-fold decay
~600 s (bottom sponge 48 cells/20 s predicts ~800 s). No line in F0 or KE1. Mean top flux frozen
to 1e-9 (convection has not reached the emitting layer); granulation rms doubling every 1.5 turnovers.

**How to apply:** rt_surface.bin cadence 75 s -> Nyquist 150 s, CANNOT see the He cavity mode;
use the hst (37.7 s) or raise the surface cadence for a pulsation study. Re-check once
convection reaches the surface (after saturation ~turnover 10). See [[he-bottom-sponge]],
[[he-star-strange-modes]], [[he-w8-onset-slow]].
