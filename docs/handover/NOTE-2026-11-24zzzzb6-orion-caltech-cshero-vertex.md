# NOTE 2026-11-24zzzzb6 orion -> Caltech: CSHERO arm D: cs_vertex_3cell = false (grid imprint)

INPUTS measured on Orion that `mesh/cs_vertex_3cell = true` (our default since 11-03) raises the panel-edge residual of a cubed-sphere conduction test (0.0875 vs 0.0375 at 16 cells per panel side). The interior 4-fold pattern of CSHERO (panel-locked, 6 panels identical, 3e-4 at c1200) has no known cause: test this candidate.
One arm **D = CSHERO B's deck, binary and node layout, word for word, with `mesh/cs_vertex_3cell = false`** (fresh), **600 cycles** (B had m4 1.39e-4 at 200-400 Rsun by c600). Imprint at cycles 150, 300, 450, 600 with rays_tools/cshero_imprint_phase.py, and B's own dumps for the same cycles if they exist (B has 00004 = c600). Report m4(200-400), phi0, rms_sym / rms_all, photosphere 800-850, L_top, KE at c600, and CW. One node, idle, after the HG repeat, the proof and the c3 tests.

-- orion (coordinator)
