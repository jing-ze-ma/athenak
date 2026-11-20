# NOTE 2026-11-20dl orion -> caltech (USER): EMA well-balanced ring test (4 arms)

The USER wants vet_rotate_avg_ema_wb tested on production before it is turned on.
New build **B_emawb** = caltech-runs/emawb, commit a6845503c889d6e550ab718cfe606cd697ace4d7. It is v20.2 plus the key, which is read only when named; the default is bitwise.
Four arms in `docs/handover/caltech_runs_1121/arms_emawb.txt`, 1 node each (NR=4), 300 cycles, **after the current 8-arm batch, as nodes free**:
- EW_on_hg vs EW_off_hg;
- EW_on_rsg vs EW_off_rsg.
All use the v20.2 production keys (+ floor_pay) with vet_rotate_avg_ema=128 named in both arms. The wb key FATALs without EMA >= 1, so if the restart's own EMA value differs from 128, please say so in the line (grep the log's parameter dump).
In each line please give:
- CW, dt, L_top / Mdot / Lz / Etot;
- the max |v| in the deep layers, if your res.py can extract it;
- the full run.log.gz and hst.

-- orion (for the USER)
