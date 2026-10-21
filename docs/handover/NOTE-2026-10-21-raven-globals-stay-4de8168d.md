# NOTE 2026-10-21 raven -> caltech + deltaai: CORRECTION — RSG E STAYS on 4de8168d (USER). Supersedes the RSG part of cdf71b98.

USER: the global runs keep their current binary 4de8168d and their decks (hesdirk2 + conservative2 + clip, rst 5 d). Do NOT switch
RSG E to 077dda23 — not now and not at the next link. If you already switched run.cfg, revert it to your 4de8168d build
(Caltech 9a56b4a1, DeltaAI c56693cb) and confirm in a NOTE.
The ACCRETOR part of cdf71b98 stands: Caltech builds 077dda23 (H200 nofma, ry_per_accretor) and starts the accretor FRESH with the
v12lzc_be_gm template (5de0635e). The accretor is the only run on 077dda23.
