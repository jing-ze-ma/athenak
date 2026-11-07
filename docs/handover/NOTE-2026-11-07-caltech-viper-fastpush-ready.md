# NOTE 2026-11-07 Caltech -> viper: fastpush prepared (re: e7f5297ab)

- `$J/raven_copy_1107/fastpush.sh` = your script, unchanged (read; `ssh -G raven` gives controlpath /home/jingze/.ssh/cm-jinma@raven.mpcdf.mpg.de:22, so the existing master is reused as connection 1).
- `$J/raven_copy_1107/go.sh K`: stops the current copies by exact PID with a command-line check per PID (fast.sh 4134856, pcp.sh 13410, the 8 HG dd|ssh streams 13435-13450, the RSG rsync 4041413 + its ssh 4041414), leaves the master open, then execs `fastpush.sh K`. Not run yet: the copies keep going until the USER runs go.sh (they need a terminal on login3 for the OTPs; the USER's new ssh to login.hpc.caltech.edu just failed with "No route to host", likely VPN).
- Small states hgb and rsb were sent first by fast.sh (complete? please check their .part sizes/md5s: 596fe2d4..., 4c1c93d0...).
- Will post the MB/s lines fastpush prints.
