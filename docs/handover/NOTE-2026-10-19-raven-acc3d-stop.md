# NOTE 2026-10-19 raven -> caltech: STOP the 3-D accretor production (v10) now (USER)

USER (10-19 ~00:15 CEST) approved: stop the Caltech acc3d production. Please scancel link 1 (4977641, running) and links 2-6
(4977642-46). KEEP the run directory, rst and dumps (do not delete). Reply with a one-line NOTE (t reached, GPU-h used).
Why: v10 over-irradiates the gainer (held donor-lobe cells emit ~cE/2 = 2x sigma T^4, and stream_sector_light = surface holds donor
light in the whole back sector): Ltop - Lmeas is 3.5x the analytic donor budget in the static phase and ~9x after stream on.
The fix (stream_sector_light = neck + donor_irr_scale ~0.35, branch acc3d-held-1018) exposes a stream-onset negative-M1-E crash
that needs the lp redesign first. A new TASK follows when the stream onset is fixed. Raven keeps a static-only copy meanwhile.
