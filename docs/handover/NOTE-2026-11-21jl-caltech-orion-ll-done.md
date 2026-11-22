# NOTE 2026-11-21jl Caltech -> orion: ll done — t9l62 head relaunched at cfl 0.6 from 00045 (USER confirmed "go" in our session)

- Head cfl-0.9 run stopped ~23:2x PST via the STOP file (graceful; it wrote t9l62.00046, moved to t9l62/rst_head09/). At stop: c~33400, CW 10, dt_retry 4 (one "last change inf").
- Relaunched the same head binary (md5 2144006b) from t9l62.00045 with the jk keys except time/cfl_number=0.6 (conv 1e-5, flim false, fwb/fcone/fwb_stage true, output5/dt 21600, VGD_XFE=0). END 00:20 PST.
- Per ll item 4: if it falls back in c33370-33460 or shows "last change inf", we stop and revert to B_j3 from 00045 (USER rule). First cycles follow.

-- Caltech
