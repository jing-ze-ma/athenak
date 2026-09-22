#!/bin/bash -l
# rstgate.sh <tag> : the postmerge_rst.sh restart gate (box mode 3, production warm
# start) run with build_clean_box_convection.  50 cycles straight vs 25 + restart + 25.
module purge >/dev/null 2>&1; module load gcc/14 cmake/4.0 >/dev/null 2>&1
A=/viper/u2/jinma/ATHENAK/athenak; G=$A/tests_cleanup_0922
X=$A/build_clean_box_convection/src/athena
IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
O="mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 meshblock/nx2=8 meshblock/nx3=8"
D=$G/rst_$1; rm -rf $D; mkdir -p $D/s $D/c
cd $D/s; timeout 1800 $X -i $IN $O time/nlim=50 > run.log 2>&1; echo "straight exit=$?"
cd $D/c; timeout 1800 $X -i $IN $O time/nlim=25 output5/dt=1e9 > run1.log 2>&1
R=$(ls rst/*.rst | tail -1); echo "restart from $R"
timeout 1800 $X -r $R time/nlim=50 > run2.log 2>&1; echo "restart exit=$?"
python3 - <<'PY'
for f in ['feczrt.hydro.hst','feczrt.user.hst']:
    a=[l for l in open('../s/'+f) if not l.startswith('#')]
    b=[l for l in open(f) if not l.startswith('#')]
    sb=set(b); miss=[i for i,l in enumerate(a) if l not in sb]
    print(f,'straight rows',len(a),'chain rows',len(b),
          'straight rows missing bit-for-bit in chain:',len(miss), miss[:5])
PY
for f in rt_surface.bin rt_profile.bin; do
  cmp -s $f ../s/$f && echo "$f IDENTICAL" || echo "$f DIFFER"
done
rm -rf $D/s/rst $D/c/rst
echo RST_GATE_DONE $1
