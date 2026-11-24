# NOTE 2026-11-24zzzd viper -> Caltech (+ orion): USER: RSG rc3 did not start from the latest restart; rerun from rsg.01809

The USER (11-24): Caltech's rc3 RSG continuation started from **rsg.01707**, which is not the latest RSG restart. **Rerun it from
the latest one, rsg.01809.**

- The latest restart is the old-scheme RSG on viper (oldprod1121, prod-1109 binary 864d120d4), which continued past rsg.01707:
  `/viper/ptmp2/jinma/orion_stag2/oldprod1121/rsg/rst/rsg.01809.rst`
  - 46402993637 B, md5 `b821975950ea1a8bb6d31009cd398d3a` (also in `rsg.01809.rst.md5` next to it)
  - header: t = 7.65731043e8 s = 8862.6 d, cycle 817373, tlim 1e30; final rst of run_04 (wall-clock stop, rc=0, 11-24 04:21 PST)
  - output last_time in the file: output1 765731000, output2/3 765676800, output4 765700000, output5/6 765504000
  - The old run's keys for that segment: `.../oldprod1121/rsg/run_04.log.cmd` (and `rsg.keys`).
- Copy it by the USER's ssh route, as you did rsg.01707, and check the md5 on arrival.
- Restart the rc3 RSG from it with the same rc3 deck/keys as the 01707 continuation (same L_bot comoving note; check the
  output last_time values against the new start). Post the job id, first cycles and s/cyc.
- The current rc3 run from rsg.01707 is superseded. Stop it once the 01809 run is up, unless the USER says otherwise.
- The RSG page build (zzzb) should use the 01809 continuation. Its history: v14 (to 8361 d), then the old-scheme segment
  01707 -> 01809 (8361 -> 8863 d, viper run_01..04 logs + `rsg.hydro.hst`/`rsg.user.hst` in oldprod1121/rsg), then rc3.

-- viper
