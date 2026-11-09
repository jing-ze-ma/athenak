#!/bin/bash
# lintol-1109 Orion CPU runs (TASK-2026-11-09-orion-lintol.md). usage: run_lintol.sh BIN OUTDIR [NR=4]
# Fresh R136 small deck (sp, 831x8x8, 4 MeshBlocks), all arms to tlim 120 s (8 cycles), keys REPLACED in deck copies.
set -e
BIN=$1; G=$2; NR=${3:-4}
H=$(cd "$(dirname "$0")" && pwd)            # docs/handover/lintol-1109
mkdir -p $G/in && cd $G || exit 1
[ -f $G/in/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt ] || gunzip -c $H/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt.gz > $G/in/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt
[ -f $G/in/planck_tops_x0.58_z0.0064_a09_T7.txt ] || gunzip -c $H/planck_tops_x0.58_z0.0064_a09_T7.txt.gz > $G/in/planck_tops_x0.58_z0.0064_a09_T7.txt
cp $H/ic_r136a1_B_dfl17.txt $G/in/
sed "s#LINTOL_DIR#$G/in#" $H/r136_small.athinput > $G/in/base.athinput
mk() {   # mk NAME key=value ... : deck copy, existing key lines replaced, new keys added to <rad_m1>
  python3 - $G/in/base.athinput $G/in/$1.athinput "${@:2}" <<'PY'
import sys, re
inp, out, kv = sys.argv[1], sys.argv[2], sys.argv[3:]
s = open(inp).read(); add = []
for x in kv:
    k, v = x.split("=", 1)
    s, n = re.subn(r"(?m)^ *%s *=.*$" % re.escape(k), "%s = %s   # lintol override" % (k, v), s)
    if n == 0: add.append("%s = %s" % (k, v))
if add:
    s, n = re.subn(r"(?m)^<rad_m1>[ \t]*$", "<rad_m1>\n" + "\n".join(add), s, count=1); assert n == 1
open(out, "w").write(s)
PY
}
SK="mesh/nx2=8 mesh/nx3=8 meshblock/nx2=4 meshblock/nx3=4 output1/dt=1e12 output2/dt=1e12 output3/dt=1e12 output4/dt=1e12 output5/dt=1e12 output6/dt=1e12 time/tlim=120.0 time/nlim=400"
run() { local n=$1; mkdir -p $G/$n; [ -f $G/$n/DONE ] && return 0
  mpirun -np $NR $BIN -i $G/in/$n.athinput $SK -d $G/$n > $G/$n/run.log 2>&1 && touch $G/$n/DONE
  echo "$n rc=$? $(grep -a 'bicgstab: outer' $G/$n/run.log | cut -c1-120)"; }
mk D
mk R2 implicit_lin_ew_max=0 implicit_lin_cnorm=1e-10 implicit_tol_rel=0.001 implicit_tol_rel_max=1e-6
mk R3 implicit_lin_ew_max=0 implicit_lin_cnorm=1e-10 implicit_tol_rel=0.0001 implicit_tol_rel_max=1e-7
mk Kall implicit_lin_tolcap=1.0
mk CE0 implicit_conv_est=0
mk CEg3 implicit_conv_est_guard=3
mk CEg10 implicit_conv_est_guard=10
for a in D R2 R3 Kall CE0 CEg3 CEg10; do run $a; done
python3 $H/acc.py $G R2 R3 D Kall CE0 CEg3 CEg10 | tee $G/ACC.txt
