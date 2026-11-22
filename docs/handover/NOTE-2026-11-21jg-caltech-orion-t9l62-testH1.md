# NOTE 2026-11-21jg Caltech -> orion: USER — t9l62 continues on Caltech testH1 first; Raven LATER from a newer rst (re lf)

The Caltech -> Raven copy runs at only ~4-7 MB/s through the MPCDF gateway (28 GB = ~90 min; Raven has no route back to Caltech). Our last holder testH1 (hpc-sm-01-11, 4 H200) runs until 00:34 PST, so the USER chose: run t9l62 HERE now, then send its newer final rst to Raven.

- t9l62 relaunched 17:53 PST on testH1 from t9l62.00037 (cycle 30085), B_p1 md5 d38d411c, kp keys + gas switch 4/6, VGD_XFE=0, cfl 0.6; wall END 00:20 PST -> graceful stop writes the final rst ~00:10.
- The partial copy of t9l62.00037 to raven:/ptmp/jinma/t9l62_from_caltech/ was STOPPED (incomplete, ignore it; nothing was written into your src_rst path). Please do NOT relaunch t9l62 on Raven from 00037.
- After the final rst we push it to raven:/raven/ptmp/jinma/stag2/prod1121/t9l62/src_rst/<name>.part, rename when complete, and post the md5 (~90 min unless the USER sets up Globus).

-- Caltech
