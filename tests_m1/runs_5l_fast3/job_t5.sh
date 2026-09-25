#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/fast3/runs/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/fast3/runs/log.err.%j
#SBATCH -J fast3_t5
#SBATCH -p apudev
#SBATCH --constraint="apu"
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# m1-fast3 timing, 1 GPU: base vs n3 (default), op_team_red, +krylov_dev, eos check_every 10
source /viper/ptmp2/jinma/fast3/scripts/common.sh
md5sum $W/bin/athena_*
N="time/nlim=60"
IN=$W/inp
for r in 1 2 3; do
  run t5 b_b0_$r 1 0 base box_convection $IN/box_vsc.athinput $N
  run t5 b_n0_$r 1 0 n3 box_convection $IN/box_vsc.athinput $N
  run t5 b_nt_$r 1 0 n3 box_convection $IN/box_vsc_t.athinput $N
  run t5 b_ntk_$r 1 0 n3 box_convection $IN/box_vsc_tk.athinput $N
  run t5 b_ne_$r 1 0 n3 box_convection $IN/box_vsc_e.athinput $N
  run t5 b_ntke_$r 1 0 n3 box_convection $IN/box_vsc_tke.athinput $N
  run t5 w_b0_$r 1 0 base none $IN/hewedge_o.athinput $N
  run t5 w_n0_$r 1 0 n3 none $IN/hewedge_o.athinput $N
  run t5 w_nt_$r 1 0 n3 none $IN/hewedge_o_t.athinput $N
  run t5 w_ntk_$r 1 0 n3 none $IN/hewedge_o_tk.athinput $N
done
run t5 pb_nt 1 1 n3 box_convection $IN/box_vsc_t.athinput time/nlim=40
run t5 pw_nt 1 1 n3 none $IN/hewedge_o_t.athinput time/nlim=40
echo ALL DONE
