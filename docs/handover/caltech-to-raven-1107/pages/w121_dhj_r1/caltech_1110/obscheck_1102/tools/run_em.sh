#!/bin/bash
# all em3d passes in one GPU session: W121 + W76 dumps 596-600 (rot 298-300), + validation (iso, fine grid)
O=/resnick/groups/carnegie_poc/jingze/obscheck_1102; cd $O
EOS=/resnick/groups/carnegie_poc/jingze/w121r1_1101/anchor/eos_1x.txt
PY=$O/venv/bin/python; T=tables/opac_R200.npz
S0=$(date +%s)
for pl in w121:/resnick/groups/carnegie_poc/jingze/w121r1_1101/runs/p1 w76:/resnick/groups/carnegie_poc/jingze/r1planets_1101/wasp76b/runs/p1; do
  k=${pl%%:*}; R=${pl#*:}
  for n in 00596 00597 00598 00599 00600; do
    [ -s out/em_${k}_$n.npz ] || $PY tools/em3d.py --eos $EOS --table $T --out out $R/bin/dhj.hydro_w.$n.bin tag=${k}_$n
  done
done
D=/resnick/groups/carnegie_poc/jingze/w121r1_1101/runs/p1/bin/dhj.hydro_w.00600.bin
[ -s out/em_val_iso.npz ] || $PY tools/em3d.py --eos $EOS --table $T --out out $D tag=val_iso nph=4 iso=2500
[ -s out/em_val_fine.npz ] || $PY tools/em3d.py --eos $EOS --table $T --out out $D tag=val_fine nb=192 nth=128 nsub=6
# cloud-deck post-processing grid (tables/cloud_cases.json), 20 phases, same dumps
for pl in w121:/resnick/groups/carnegie_poc/jingze/w121r1_1101/runs/p1 w76:/resnick/groups/carnegie_poc/jingze/r1planets_1101/wasp76b/runs/p1; do
  k=${pl%%:*}; R=${pl#*:}
  for n in 00596 00597 00598 00599 00600; do
    [ -s out/emcl_${k}_$n.npz ] || $PY tools/em3d.py --eos $EOS --table $T --out out $R/bin/dhj.hydro_w.$n.bin tag=${k}_$n clouds=tables/cloud_cases.json nph=20
  done
done
echo "ALL em3d wall $(( $(date +%s)-S0 )) s"
