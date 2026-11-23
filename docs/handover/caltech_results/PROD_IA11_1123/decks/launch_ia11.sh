#!/bin/bash
# launch_ia11.sh STAR NODE HOLDER END   (orion NOTE 2026-11-23i, USER "stop IA9 and launch IA11": FRESH bsg/m35mms/t9l62 on IA11 = caltech-ia11-bd9553537, decks caltech_ia11 with Caltech file paths; NO keys on the command line)
# 1 node = 4 H200 (NR=4), as the 11-20/21 productions of each run. Fresh start while $O/$STAR/rst has no .rst; same command later = restart from the newest rst.
# END = "YYYY-MM-DD HH:MM" (holder end - 20 min). Env modelled on switch_p1.sh (saved env first, THEN NODE/HOLDER/OUT).
STAR=$1; N_=$2; H_=$3; E_=$4; : "${STAR:?}" "${N_:?}" "${H_:?}" "${E_:?}"
P=/resnick/groups/carnegie_poc/jingze/prod_fresh_1121b; O=/resnick/groups/carnegie_poc/jingze/prod_ia11_1123
while IFS= read -r l; do export "$l" 2>/dev/null; done < $P/i1_m35mms_f1k.env
export NODE=$N_ HOLDER=$H_ END="$E_" OUT=$O KEEPRST=3 FRESH=/resnick/groups/carnegie_poc/jingze/fresh_prod_1120 MEMWRAP=/resnick/groups/carnegie_poc/jingze/orion_runs_1121/ctl/memwrap.sh
export VGD_XFE=0   # as switch_p1.sh (also the code default since rsthip-1121)
BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia11_bd9553537
exec /resnick/groups/carnegie_poc/jingze/fresh_prod_1120/fg_v2/prodg.sh $STAR $BIN $O/decks_final/IA11_$STAR.athinput
