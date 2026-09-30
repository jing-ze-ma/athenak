#!/bin/bash
# all configurations on the login node, <= 16 processes, ~0.4 GB each
cd /viper/ptmp2/jinma/w121prod_0929/synth2_trans
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R="nice -n 10 python3 trans.py run"
$R 10x template wind=0 rot=0 ep=mid &
$R 1x base
wait
$R 10x base
$R 1x ext ext=1
$R 10x ext ext=1
$R 1x s5base nep=5 & $R 1x s5fe01 nep=5 fe=0.1 & $R 1x s5fe10 nep=5 fe=10 & wait
$R 1x s5ext nep=5 ext=1 & $R 1x s5pext01 nep=5 ext=1 pext=0.1 & $R 1x s5pext10 nep=5 ext=1 pext=10 & wait
$R 1x s5norot nep=5 rot=0 & $R 1x s5nowind nep=5 wind=0 & $R 1x s5gas nep=5 chem=gas & wait
echo ALLDONE
