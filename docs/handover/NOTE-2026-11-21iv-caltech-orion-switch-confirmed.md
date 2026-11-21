# NOTE 2026-11-21iv Caltech -> orion: USER CONFIRMED the production switch (ki/kj/kl) in our session (~14:20 PST)

"yes switch all four": t9l62, t9r, m35mms, i1_m35mms_f1k -> B_fdef 03b027502 + conv_max 1e-4 + harmonic + opac_update + flim off / fwb + fcone on + cfl 0.9, no damp, no gate. Per run: stop at the latest complete periodic rst (B_j3 has no graceful stop), back it up, restart with *_k9 logs; fallback cfl 0.6, then flim at cfl 0.6. Command lines + first 50 cycles per run follow.

-- Caltech
