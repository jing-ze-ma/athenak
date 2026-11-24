# NOTE 2026-11-24zzzs orion -> Caltech: thinf pairs on the cs stars too (USER: set it true everywhere if harmless)

The USER wants `clean_stag_thinf=true` on all stars if it makes no difference. The sp stars: M35ZAMS, m35mms, bsg harmless (100-300 cycles); t9l62 pending. Add the cs ones when nodes allow, as cheap pairs (true vs false, same rst, key on the command line only, 100 cycles, hst + log, every key named):
1. **RSG** from rsg.01707 or, once its md5 is checked, rsg.01809 (rc3 binary and its keys; one pair; 3 nodes, so after HG and the t9l62 pair).
2. **HG** from hegiant.00205: run H0 with thinf=true as one extra 51-cycle arm (HT) beside H0 / H1 of 24zzzr; only if a node is free, after H0 and H1.
Report dt, L_top / L_mid / L_int / L_bot, passes, its, CW / fallback / retry, s/cyc, for each pair.
R136: no H200 binary/rst arm is needed now.

-- orion (coordinator)
