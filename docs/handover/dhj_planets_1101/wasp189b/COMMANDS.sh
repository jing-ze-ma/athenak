#!/bin/bash
# WASP-189b package: the command lines that produced every file (run from the package dir, login node, nice).
# Not meant to be re-run blindly: each step was checked before the next.
W=/viper/ptmp2/jinma/dhj_planets_1101/wasp189b; PY=/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2025.06/bin/python3
cd $W || exit 1
# SED (ATLAS9 ckp02/ckp05 8000/8250 K from archive.stsci.edu cdbs/grid/ck04models into sed/atlas9)
$PY sed/make_sed.py 8000 4.00 sed/sw_band_flux_W189_11.txt 0.29        # -> sed/make_sed.out
$PY sed/make_sed.py 8000 4.00 sed/sw_band_flux_W189_11_p00.txt 0.0     # [M/H] 0 comparison
# 3x: face-0 RCE (guess G0/AP) -> face-0 at planet_setup's ap/grav -> f5 (RESTOL 2e-5) -> transit, IC, setup
bash chain.sh $W/rce/rceA 2122 1.09e10 /viper/ptmp2/jinma/w121prod_0929/met3/rceB3f/ckrce_nq2_B.npz 400 nof5
bash chain.sh $W/rce/rceB 2109.824 1.093122e10 $W/rce/rceA/ckrce_nq2_B.npz 400     # face 0 done; f5 stalled at res ~4e-6
cp rce/rceB/f5/ckpt.npz rce/rceB/f5/ckpt_it40.npz
bash finish_f5.sh $W/rce/rceB 2109.824 1.093122e10 $W/rce/rceB/f5/ckpt_it40.npz
cp rce/rceB/ic.txt ic/ic_w189_3x.txt
# dayside 1-D estimates (table coverage only)
scripts/run_rce_day.sh $W/rce/day_mean 2109.824 1.093122e10 $W/rce/rceA/ckrce_nq2_B.npz 3136.6
scripts/run_rce_day.sh $W/rce/day_sub 2109.824 1.093122e10 $W/rce/rceA/ckrce_nq2_B.npz 3730.1
# 1x
W189_MET=1x bash chain.sh $W/rce/rce1A 2109.824 1.093122e10 $W/rce/rceA/ckrce_nq2_B.npz 400 nof5
#   rce1B face 0 (chain.sh) stopped by me at it 139 (res 1.7e-6); its it-120 checkpoint seeded f5:
cp rce/rce1B/ckpt_nq2_B.npz rce/rce1B/ckpt_face0_it120.npz
W189_MET=1x bash finish_f5.sh $W/rce/rce1B 2092.501 1.097637e10 $W/rce/rce1B/ckpt_face0_it120.npz
cp rce/rce1B/ic.txt ic/ic_w189_1x.txt
# grids (W121 hfresh cells/H carried over), cost
(cd grid/g3x && $PY $W/scripts/grid112.py $W/ic/ic_w189_3x.txt 2118.7772 1.090810e10 112 120 124 128 > grid.log)
(cd grid/g1x && G112_EOSF=/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt \
   $PY $W/scripts/grid112.py $W/ic/ic_w189_1x.txt 2100.830 1.095459e10 124 > grid.log)
$PY scripts/cost.py ic/ic_w189_3x.txt 2118.7772 1.090810e10 grid/g3x/faces_124.txt 235356.5 W189_3x_n124 > gate/cost_3x_n124.txt
$PY scripts/cost.py ic/ic_w189_3x.txt 2118.7772 1.090810e10 grid/g3x/faces_112.txt 235356.5 W189_3x_n112 > gate/cost_3x_n112.txt
G112_EOSF=/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt \
  $PY scripts/cost.py ic/ic_w189_1x.txt 2100.830 1.095459e10 grid/g1x/faces_124.txt 235356.5 W189_1x_n124 > gate/cost_1x_n124.txt
# decks
python3 scripts/build_deck.py 3x rce/rceB/setup.txt grid/g3x/gridkeys_124.txt deck/w189_3x.athinput
python3 scripts/build_deck.py 1x rce/rce1B/setup.txt grid/g1x/gridkeys_124.txt deck/w189_1x.athinput
# table coverage
$PY scripts/table_coverage.py > gate/table_coverage_tables.txt
TC_MET=1x $PY scripts/table_coverage.py > gate/table_coverage_tables_1x.txt
TC_MET=10x $PY scripts/table_coverage.py > gate/table_coverage_tables_10x.txt
$PY scripts/table_coverage.py ic_3x=ic/ic_w189_3x.txt global_f5_3x=rce/rceB/f5/ckrce_nq2_B.txt \
  dayside_mean_3x=rce/day_mean/ckrce_nq2_B.txt substellar4x_3x=rce/day_sub/ckrce_nq2_B.txt > gate/table_coverage_3x.txt
TC_MET=1x $PY scripts/table_coverage.py ic_1x=ic/ic_w189_1x.txt global_f5_1x=rce/rce1B/f5/ckrce_nq2_B.txt > gate/table_coverage_1x.txt
# CPU gates (binary /viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_8bcd50c8)
(cd gate && bash gate.sh 3x && bash gate.sh 1x)
(cd gate && $PY gate_ic.py run_gz_3x/bin/dhj.hydro_w.00000.bin ../ic/ic_w189_3x.txt ../grid/g3x/gridkeys_124.txt \
   1.195e-3 1.156186e10 1.158728e10 2118.7772 2.669647e-05 > gate_ic_3x.txt)
(cd gate && GATE_EOSDIR=/viper/ptmp2/jinma/deepconv_0925/dump $PY gate_ic.py run_gz_1x/bin/dhj.hydro_w.00000.bin \
   ../ic/ic_w189_1x.txt ../grid/g1x/gridkeys_124.txt 2.369e-3 1.156186e10 1.158728e10 2100.8301 2.669647e-05 > gate_ic_1x.txt)
# Roche anchors (fork/fin51-r-fix c8aac483c roche_anchor.py), dynamo, system/tide
(cd gate && $PY ../scripts/roche_anchor.py --mp 1.99 --rp 1.619 --pref 1.195e-3 --a_au 0.050529 --porb 2.7240330 \
   --teff 8000 --rstar 2.362 --mstar 2.030 --profile ../ic/ic_w189_3x.txt \
   --eos /viper/ptmp2/jinma/w121prod_0929/met3/dump/eos_table.txt --x1max 1.567297e10 > roche_anchor_3x.txt)
(cd gate && $PY ../scripts/roche_anchor.py --mp 1.99 --rp 1.619 --pref 2.369e-3 --a_au 0.050529 --porb 2.7240330 \
   --teff 8000 --rstar 2.362 --mstar 2.030 --profile ../ic/ic_w189_1x.txt \
   --eos /viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt --x1max 1.595318e10 > roche_anchor_1x.txt)
python3 scripts/dynamo.py > gate/dynamo.txt
python3 scripts/system.py 1.090810e10 1.157455e10 1.567297e10
# GPU smoke (apudev, hfresh binary, both decks 40 cycles): sbatch smoke/smoke.sh
# Production at Caltech (run dir with ck -> package ck/, ic -> package ic/):
#   athena -i deck/w189_1x.athinput -t <wall>    (nominal)   |   deck/w189_3x.athinput (sensitivity / trend-consistent)
