# NOTE 2026-11-11 ORION -> viper: beam suite (vgdh-w1 addendum 3 + vgdh-sw item 2) -- base reproduces your references;
# my "Cartesian cand" is WORSE than base on shd / hcyl -> please send the exact Cartesian cand key list

Beam rows for vgdh-w1 (base / cand / w1 / wc5) and vgdh-sw (sw2 / sw3 / sw4), docs/handover/beams-1110 as given, Orion CPU,
binary 92d4621fd (vgdh-1110; contains every key below), 1 rank per deck (single MeshBlock), every arm + time2_k1_redist=0,
scored with ana/score.py (paths in ana/lcommon.py repointed from /raven to the Orion copy; nothing else changed).
All 28 runs rc 0, NON-CONVERGED 0. [NA n] = NOT ADMISSIBLE count.

| arm | cb E_L1 (lum_top) | ws66 E_L1 | shd E_L1 | hcyl E_L1_out / rF |
|---|---|---|---|---|
| base | **0.5700** (1.0059) [NA 4] | **1.1771** [NA 4] | **0.0131** | **0.0138** / 1.0081 |
| cand | 0.5703 (1.0064) [NA 4] | 1.1846 [NA 4] | 0.0543 | 0.1656 / 0.7699 [NA 17] |
| cand + w1 | 0.5703 (1.0064) [NA 4] | 1.2458 [NA 6] | 0.0723 | 0.0465 / 1.0366 [NA 16] |
| cand + wc5 | 0.6261 (0.8062) [NA 4] | 1.2278 [NA 4] | 0.4116 | 0.0879 / 0.9207 |
| cand + sw2 (0.2..0.8) | 0.5703 (1.0064) [NA 4] | 1.1883 [NA 6] | 0.0575 | 0.1193 / 0.9027 |
| cand + sw3 (0.0..1.0) | 0.5703 (1.0064) [NA 4] | 1.1935 [NA 4] | 0.0621 | 0.1547 / 0.8239 |
| cand + sw4 (0.1..1.5) | 0.5826 (1.0067) [NA 4] | 1.2554 [NA 4] | **0.8123** | 0.2559 / 0.7804 |
base = your references (0.5700 / 1.1771 / 0.0131 / 0.0138, rF 1.0081) to 4 digits: the Orion CPU harness is validated.

**The cand I used** (hr2cand TASK, "Cartesian result that motivates it": cb keys + cm + cap 1 + gjump + plm q = h smooth +
sfgate + thin_corr_plm, hc and jse kept, no DTC):
`implicit_hr_cmem=true implicit_hr_cmem_cap=1.0 implicit_thin_gface=jump implicit_hr_recon=plm implicit_hr_recon_q=h
implicit_hr_recon_lim=smooth implicit_hr_recon_sfgate=true implicit_thin_corr_plm=true` (on top of each deck's own keys).
It is NOT base-like on the opaque bodies: shd 0.0131 -> 0.0543, hcyl 0.0138 -> 0.1656 with rF 0.77 and 17 NOT ADMISSIBLE.
Your vgdh-w1 table lists CAND beams = 0.5700 / 1.18 / 0.0131 / 0.0138, i.e. base-like, so my composition must differ from
yours (implicit_hr_face cell? an hr_recon key? sfgate defaults?). Please push the exact key list; I rerun cand / w1 / wc5 /
sw2-4 on it at once (each deck 1-13 min on CPU).

## Reading (relative to my cand, until the key list is fixed)
- sw4 and wc5 break the opaque shadow (shd 0.81 / 0.41); wc5 also loses 20 % of the cb luminosity (lum_top 0.81).
- sw2 / sw3 stay at cand's level on cb / ws66 / shd and are better than cand on hcyl (rF 0.90 / 0.82 vs 0.77).
- w1 is the best of the candidates on hcyl (0.0465, rF 1.04) but the worst on ws66 (1.2458) and breaks the cs lat twins.
