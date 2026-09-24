#!/bin/bash
# final T-S2 and T-S5 runs with the final binary
cd /viper/ptmp2/jinma/s2_0924
./run.sh F_fs32 new none 1 $PWD/inp/sph_fs.athinput mesh/nx1=32 meshblock/nx1=32 &
./run.sh F_fss32 new none 1 $PWD/inp/sph_fs_str.athinput mesh/nx1=32 meshblock/nx1=32 &
./run.sh F_fs64 new none 1 $PWD/inp/sph_fs.athinput mesh/nx1=64 meshblock/nx1=64 &
./run.sh F_fss64 new none 1 $PWD/inp/sph_fs_str.athinput mesh/nx1=64 meshblock/nx1=64 &
./run.sh F_fs128 new none 1 $PWD/inp/sph_fs.athinput mesh/nx1=128 meshblock/nx1=128 &
./run.sh F_fss128 new none 1 $PWD/inp/sph_fs_str.athinput mesh/nx1=128 meshblock/nx1=128 &
./run.sh F_rwc10_none new none 1 $PWD/inp/rw_cart.athinput rad_m1/implicit_offdiag=none mesh/x2min=0.20796326794895137 mesh/x2max=1.2079632679489514 &
./run.sh F_rws10_none new none 1 $PWD/inp/rw_sph_R10.athinput rad_m1/implicit_offdiag=none &
./run.sh F_rwc100_none new none 1 $PWD/inp/rw_cart.athinput rad_m1/implicit_offdiag=none mesh/x2min=0.5796326794930167 mesh/x2max=1.5796326794930167 &
./run.sh F_rws100_none new none 1 $PWD/inp/rw_sph_R100.athinput rad_m1/implicit_offdiag=none &
./run.sh F_rwc100_lagged new none 1 $PWD/inp/rw_cart.athinput rad_m1/implicit_offdiag=lagged mesh/x2min=0.5796326794930167 mesh/x2max=1.5796326794930167 &
./run.sh F_rws100_lagged new none 1 $PWD/inp/rw_sph_R100.athinput rad_m1/implicit_offdiag=lagged &
./run.sh F_rwc1000_none new none 1 $PWD/inp/rw_cart.athinput rad_m1/implicit_offdiag=none mesh/x2min=0.2963267950694899 mesh/x2max=1.29632679506949 &
./run.sh F_rws1000_none new none 1 $PWD/inp/rw_sph_R1000.athinput rad_m1/implicit_offdiag=none &
./run.sh F_rwc10000_none new none 1 $PWD/inp/rw_cart.athinput rad_m1/implicit_offdiag=none mesh/x2min=0.4632679158166866 mesh/x2max=1.4632679158166866 &
./run.sh F_rws10000_none new none 1 $PWD/inp/rw_sph_R10000.athinput rad_m1/implicit_offdiag=none &
wait
