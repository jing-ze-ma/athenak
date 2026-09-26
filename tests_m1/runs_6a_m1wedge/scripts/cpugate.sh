#!/bin/bash
# CPU correctness gates for sph_wedge: restart bitwise (He + grey), 2 vs 4 ranks (He),
# a 300-cycle He relaxation run, T-S4 bitwise base vs new.  usage: cpugate.sh <new tag>
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/sprhd_0926; B=$W/bin/athena_$1_none_cpu; C=$W/cpu/gate_$1; mkdir -p $C
r() { local d=$1 np=$2; shift 2; mkdir -p $C/$d; (cd $C/$d && timeout 3000 mpirun --oversubscribe -np $np "$@" -d . > log.txt 2>&1; echo "rc=$?" >> log.txt); }
HE="$W/inp/he.athinput"; GA="$W/inp/gateA.athinput"
( r he_full 4 $B -i $HE time/nlim=40 output4/dcycle=20
  r he_rr 4 $B -r $C/he_full/rst/hw.00001.rst time/nlim=40 ) &
( r he_np2 2 $B -i $HE time/nlim=40 output4/dcycle=20 ) &
( r ga_full 4 $B -i $GA time/nlim=40 output4/dcycle=20 mesh/nx2=8 mesh/nx3=8 meshblock/nx2=4 meshblock/nx3=4
  r ga_rr 4 $B -r $C/ga_full/rst/wg.00001.rst time/nlim=40 ) &
( r he_long 4 $B -i $HE time/nlim=300 output2/dt=10.0 output3/dt=10.0 ) &
TS=/viper/ptmp2/jinma/vetevery_0926/wedge.athinput
TG="mesh/nx2=32 mesh/nx3=32 meshblock/nx2=16 meshblock/nx3=16 time/nlim=12"
( r ts4_base 4 $W/bin/athena_base_none_cpu -i $TS $TG
  r ts4_new 4 $B -i $TS $TG ) &
wait
python3 - <<PY
import glob
def body(f):
    a=open(f,'rb').read(); return a[a.find(b'<par_end>'):]
C='$C'
for p,b in [('he','hw'),('ga','wg')]:
    fa=sorted(glob.glob(f'{C}/{p}_full/rst/{b}.*.rst'))[-1]; fb=sorted(glob.glob(f'{C}/{p}_rr/rst/{b}.*.rst'))[-1]
    print(p,'restart bitwise:', body(fa)==body(fb), fa.split('/')[-1], fb.split('/')[-1])
import filecmp,os
for f in sorted(os.listdir(f'{C}/ts4_base')):
    if f.endswith('.txt'): continue
    pa=f'{C}/ts4_base/{f}'; pb=f'{C}/ts4_new/{f}'
    if os.path.isdir(pa):
        for g in sorted(os.listdir(pa)): print('T-S4',f+'/'+g, filecmp.cmp(pa+'/'+g,pb+'/'+g,shallow=False))
    else: print('T-S4',f, filecmp.cmp(pa,pb,shallow=False))
PY
echo CPUGATE DONE
