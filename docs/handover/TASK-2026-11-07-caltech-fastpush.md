# TASK 2026-11-07 viper -> Caltech (USER: "can you just set it for me"): one-command fast transfer; please prepare it for the USER

Viper cannot reach Caltech (no VPN), and every new MPCDF connection needs the USER's OTPs, so the USER runs ONE command; you prepare it.
Measured (NOTE-2026-11-07-viper-caltech-transfer-separate-connections, b75707e95): everything shares ONE connection, ~9 MB/s total.

Please, on the Caltech login node that holds the current Raven master:
1. Copy `docs/handover/caltech-to-raven-1107/fastpush.sh` (this branch) to `$J/raven_copy_1107/fastpush.sh`, chmod +x. Read it first.
2. Stop the current copies **by exact PID only** (read each full command line first; never pkill patterns): the `fast.sh` HG streams (dd | ssh ... hegiant.00139.rst.part) and the single rsync of rsg.01295.rst (and its parent if any). fastpush.sh re-sends both files into the same names (RSG to `.part`; the HG byte ranges overwrite identically, so the partial HG .part is fine). Leave the Raven master connection itself OPEN (fastpush uses it as connection 1).
3. Tell the USER (one line): "in tmux on <this login node>: `bash $J/raven_copy_1107/fastpush.sh 3` — type gate OTP then raven OTP three times". If the USER prefers fewer OTPs: `fastpush.sh 1` (one extra connection) — the printed MB/s after RSG tells whether more connections help.
4. Post a NOTE with the MB/s lines fastpush prints. Viper watches both .part files (md5 vs HANDOVER sec 0) and renames them; viper's held RSG/HG smokes start after that.
