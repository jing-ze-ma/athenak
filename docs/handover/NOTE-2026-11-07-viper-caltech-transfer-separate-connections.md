# NOTE 2026-11-07 viper -> Caltech (USER: "can you make the transfer faster from caltech"): 8 streams on ONE connection give no gain; use separate TCP connections

## Measured on Raven (allocated bytes over 60 s, ~12:00 PST)
- HG `hegiant.00139.rst.part` (fast.sh, 8 streams): 12.8 of 97.5 GB written, **~6 MB/s** total.
- RSG rsync (`.rsg.01295.rst.*`, single stream): 23.2 of 45.4 GB, ~0 MB/s in that minute.
- Sum stays ~9 MB/s = the old single-rsync rate. All streams share ONE ssh master, and ProxyJump sends the master through ONE gate connection: the per-connection limit is the cap, not the stream count. At this rate the remaining ~107 GB need > 3 h.

## Fix: K independent connections (each its own gate hop). The USER opens them (gate OTP + raven OTP each); you then spread the streams.
On the Caltech login node, for k = 2 3 4 (USER, one at a time):
```
ssh -o ControlMaster=yes -o ControlPersist=12h -o ControlPath=~/.ssh/cmk$k-%r@%h:%p \
    -o ProxyCommand="ssh -o ControlMaster=no -o ControlPath=none -W %h:%p mpcdf-gate" -fN raven
```
(`ControlPath=none` on the gate hop is what makes it a NEW gate connection instead of reusing cm-...@gate1.)
Then pcp with each stream on its own connection, round robin: in pcp.sh replace `ssh -o Compression=no raven` by
`ssh -o Compression=no -o ControlPath=~/.ssh/cmk$(( i % K + 1 ))-%r@%h:%p raven` (k=1 = the existing master: use its ControlPath).
1. First a 1-minute test: 1 stream on cmk2 alone vs the existing load -> if the sum rises above ~9 MB/s, the per-connection cap is confirmed.
2. **Priority: RSG first** (both RSG smokes on Raven/viper wait for it; viper smokes are held until it lands). Re-send RSG with pcp over all K connections to `raven_gap_1105/rsgE/rst/rsg.01295.rst.part` (viper's watcher md5-checks .part files and renames them; it will be extended to RSG), then stop the single rsync of RSG. Then HG: continue the running fast.sh streams, or restart HG over all K connections into the SAME `.part` name (byte ranges overwrite identically; a restart from scratch is safe).
3. Post a NOTE with the measured rates.
If MPCDF / Caltech offer Globus endpoints for /raven/ptmp and your $J, that would beat all of this; we have not checked (USER decision).
