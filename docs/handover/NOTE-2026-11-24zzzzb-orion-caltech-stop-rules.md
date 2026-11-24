# NOTE 2026-11-24zzzzb orion -> Caltech: priorities and stop rules (USER: nothing unnecessary)

Order by what holds others up: (1) HG H0 / H1 / HT (the method + cs streaming; blocked on hegiant.00205), (2) S5, S3, S2, S1 (fresh-star settings), (3) 24zzzza items 1-3 (logs, 7x36, panel-edge number), (4) ZA / ZB / ZC, (5) CPU items zzzu, zzzo.
Stop rules, so that no GPU time goes to settled questions:
- An `on` arm that shows a FAIL of 24zzzy (extra CW / fallback / BE / retry, any ef<=0): stop it then, post it, and don't start the remaining S arms until orion answers.
- H1 hits its early-stop rule (24zzzr: dt < 45 s for 5 cycles, or the first stage fallback): same, and no more S arms.
- S4 (t9l62 quiet) runs only if S5 passes. Skip ZC if ZB and ZA both FATAL or the node is needed for HG.
- HG preempts ZA / ZB / ZC and S1 / S2 when the file lands.
Please post: the node table now, and bytes copied / ETA of hegiant.00205.

-- orion (coordinator)
