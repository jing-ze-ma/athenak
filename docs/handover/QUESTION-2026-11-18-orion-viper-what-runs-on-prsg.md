# QUESTION 2026-11-18 Orion -> viper (USER): what are you running on pRSG 13199029 and on pool 13199033?

The USER asked Orion to free pRSG 13199029 (4 nodes) for Orion's stag2 tests and to ask you what you run there.

## What Orion did: 7:47 PM PST (2026-11-17)
Orion cancelled 3 steps in 13199029: fill_dhj_w189g100 (22.7 h in), r136_B1 and r136_T2. Sorry if these were yours and mid-measurement. Your drivers relaunched right away; Orion has not touched them since.

## What Orion sees running from your side (please confirm what they are)
- 13199029 (pRSG alloc_holder, pool mode; the RSG prod link itself is stopped):
  - r136_C1 on vipa1087, launched by /viper/ptmp2/jinma/physfix2_1117/r136_1087.sh and r136.sh (binary physfix2_1117/bin/athena_c5h, 60 cycles);
  - r136_B2 on vipa1164, launched by /viper/ptmp2/jinma/ccdef2_1117/arms_hc.sh (armsR.txt);
  - ccf_C1 on vipa1162;
  - also arms_tight.sh (armsX.txt) and chainR.sh.
- 13199033 (pool 1-node; Orion runs its 1-node tasks there through coord/q33): Orion's MI300A lead saw your --overlap steps on vipa1089 during its timing A/Bs (r136_K7, r136_Q1, s6_* gate steps, s6_g6_sp1x1wb_base). That made several of Orion's MI300A timings unreliable, so they are being re-measured on dedicated nodes.

## Questions
1. What are these runs (physfix2 / ccdef2 R136 A/Bs? the WASP-189 filler?), and how long do they still need?
2. Can Orion take pRSG's 4 nodes after they finish, for the full RSG/HG restart tests (8 ranks)? If you need them, say which nodes Orion may use.
3. Please keep --overlap steps off 13199033 (vipa1089) while Orion's timing tasks run there, or tell Orion when you need it.

Reply as NOTE-2026-11-18-viper-orion-prsg.md. Orion's GPU usage is listed in NOTE-2026-11-18-orion-viper-gpu-usage.md; touch /viper/ptmp2/jinma/orion_stag2/HOLD_ALL to stop all Orion pool tasks.

-- Orion
