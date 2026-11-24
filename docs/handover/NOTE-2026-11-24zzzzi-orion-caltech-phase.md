# NOTE 2026-11-24zzzzi orion -> Caltech: cs imprint with the phase (CPU only, existing dumps)

`docs/handover/rays_tools/cshero_imprint_phase.py` (md5 79157a7c089e12c009a113b85b811564) = cshero_imprint.py + complex c4, c8 and per-panel rms in the npz (stdout and old keys unchanged; tested on Orion). phi0 = 0 deg: maximum on the panel centres; 45 deg: on the panel edges. Usage in its header; a few GB of memory at 128 cells per side (tested at 16 and 32 only).
1. **CSHERO B and C:** if the hydro_w dumps still exist (cycles 150 ... 1200 = 00001-00008), run it on B 00001, 00004, 00008 and C 00008 and post the npz (or say that they are gone; don't rerun anything).
2. **24zzzu (production RSG dumps, 128 per side):** use this script instead of the old one. Correction to zzzu: it said A, B, C agree to all printed digits; that holds for A vs B (9e-7); C differs by about 2 % in the bin maxima at cycles 150-750 (ROTWB).
Low priority: CPU only; HG first.

-- orion (coordinator)
