#!/bin/bash -l
#SBATCH -o log.out.%j
#SBATCH -e log.err.%j
#SBATCH -J RG_bisect2
#SBATCH -p apu
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:40:00
# ---------------------------------------------------------------------------------------
# RG restart non-faithfulness: ARM BISECTION, ROUND 2.  See NOTES.md.
#
# Round 1 (job 11651679) turned every physics switch off in turn and ALL TEN arms still
# DIFFERed at cycle 55801.  So the unrestored per-cycle state is upstream of every one of
# them.  Round 2 therefore (a) LOOKS AT THE GHOSTS instead of guessing, and (b) removes
# the three remaining state-carrying pieces of the update (the well-balanced background,
# the etotgrav potential offset, the position-aware PLM slopes) and the MPI seam.
#
# FOUR runs per ARM, all with the arm's overrides:
#   P  : -r start.rst (cycle 55650)  -> nlim=55800 ; writes the arm's own rst at 55800
#   E  : -r start.rst (cycle 55650)  -> nlim=55801 ; continuous, one cycle further
#   F0 : -r P/rst/<last>.rst (55800) -> nlim=55800 ; restarted, ZERO cycles
#   F  : -r P/rst/<last>.rst (55800) -> nlim=55801 ; restarted, one cycle
#
# Two comparisons, both WITH GHOST ZONES (output2/ghost_zones=true):
#   P vs F0 at 55800 : "does the restart + Driver::Initialize reproduce the ghosts?"
#                      (this is arm F of the brief, the soundness control for the method:
#                       zero cycles are integrated, so any difference is the init path)
#   E vs F  at 55801 : the phenomenon itself, now resolved by ghost class.
# Cells are split into active / radial ghost / angular face ghost / angular corner ghost /
# mixed (radial x angular).  Block is 480x16x16 with nghost=3, so the dumped arrays are
# 486x22x22: active is i in [3,482], j,k in [3,18].
#
# time/ndiag=1 makes every cycle print "cycle=N time=... dt=...", so the dt the restarted
# run takes at 55800 can be read against the continuous one (see NOTES.md, candidate 1).
# ndiag, the output dt's and time/dt_min are all bitwise inert.
# ---------------------------------------------------------------------------------------
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"

BASE=/viper/u2/jinma/ATHENAK/bench/RG_fofc_long
BW=$BASE/bitwise
START=$BW/E0/start.rst          # restart at cycle 55650, t = 1.700000e6
N0=55800                        # nlim is ABSOLUTE: the rst cycle
N1=55801                        # one cycle further
OUT="output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 time/dt_min=0.0 \
output2/ghost_zones=true time/ndiag=1"

# name | overrides | nranks
ARMS=(
  "base|            |2"
  "wb_off|hydro/wellbalance_dynamic=false hydro/wb_x1=false|2"
  "etotgrav_off|hydro/etotgrav=false|2"
  "recon_dc|hydro/reconstruct=dc|2"
  "onerank|            |1"
)

for spec in "${ARMS[@]}"; do
  IFS='|' read -r name ov NR <<< "$spec"
  echo "=============== ARM $name : nranks=$NR overrides [$ov] $(date +%T)"
  for r in P E F0 F; do
    mkdir -p "$name/$r"
    ln -f "$BASE/athena" "$name/$r/athena"
    cp -f "$BW/E0/rg.athinput" "$name/$r/rg.athinput"
  done
  ln -f "$START" "$name/P/start.rst"
  ln -f "$START" "$name/E/start.rst"
  ( cd "$name/P" && srun -n $NR ./athena -r start.rst -i rg.athinput -d . \
      time/nlim=$N0 $OUT $ov > out.txt 2>&1 )
  ( cd "$name/E" && srun -n $NR ./athena -r start.rst -i rg.athinput -d . \
      time/nlim=$N1 $OUT $ov > out.txt 2>&1 )
  R=$(ls -1 "$name"/P/rst/*.rst 2>/dev/null | tail -1)
  if [ -z "$R" ]; then echo "ARM $name: NO RST FROM P -- arm FAILED, skipping"; continue; fi
  echo "  P rst: $R"
  ln -f "$R" "$name/F0/start.rst"
  ln -f "$R" "$name/F/start.rst"
  ( cd "$name/F0" && srun -n $NR ./athena -r start.rst -i rg.athinput -d . \
      time/nlim=$N0 $OUT $ov > out.txt 2>&1 )
  ( cd "$name/F" && srun -n $NR ./athena -r start.rst -i rg.athinput -d . \
      time/nlim=$N1 $OUT $ov > out.txt 2>&1 )
  for r in P E F0 F; do
    printf "  %s/%s: %s | %s\n" "$name" "$r" \
      "$(grep -m1 'Terminating on' $name/$r/out.txt || echo 'NO TERMINATION LINE')" \
      "$(ls $name/$r/bin/*.bin 2>/dev/null | tail -1)"
  done
  echo "  --- dt at the restart cycle (7 digits only, a coarse check of candidate 1)"
  for r in P E F0 F; do
    printf "    %-3s %s\n" "$r" \
      "$(grep -E 'cycle=(55799|55800|55801) ' $name/$r/out.txt | tail -3 | tr '\n' ' ')"
  done
done

echo "=============== COMPARISON $(date +%T)"
python3 - <<'PY'
import sys, glob, os, numpy as np
sys.path.insert(0,'/viper/u2/jinma/ATHENAK/athenak/vis/python'); import bin_convert as bc

NX1, NX2, NX3 = 480, 16, 16      # <meshblock> of rg.athinput

def last(d):
    fs = sorted(glob.glob(d+'/bin/*.bin'))
    if not fs: return None, None
    return fs[-1], bc.read_binary(fs[-1])

def classify(shape):
    """shape = (nmb, nk, nj, ni) -> ghost widths and the active slices."""
    _, nk, nj, ni = shape
    gi, gj, gk = (ni-NX1)//2, (nj-NX2)//2, (nk-NX3)//2
    return gi, gj, gk

def report(tag, A, B):
    """Compare every variable; split the differing cells by ghost class."""
    same = True
    for v in A['var_names']:
        da = np.stack([np.asarray(m) for m in A['mb_data'][v]])
        db = np.stack([np.asarray(m) for m in B['mb_data'][v]])
        if da.shape != db.shape:
            print(f"  {v}: SHAPE MISMATCH {da.shape} vs {db.shape}"); same = False; continue
        neq = (da != db)
        nd = int(neq.sum())
        if nd == 0:
            print(f"  {v}: bitwise identical  (shape {da.shape})")
            continue
        same = False
        gi, gj, gk = classify(da.shape)
        idx = np.argwhere(neq)                       # columns: m, k, j, i
        _, nk, nj, ni = da.shape
        ig = (idx[:,3] < gi) | (idx[:,3] >= gi+NX1)  # radial ghost
        jg = (idx[:,2] < gj) | (idx[:,2] >= gj+NX2)
        kg = (idx[:,1] < gk) | (idx[:,1] >= gk+NX3)
        ang = jg.astype(int) + kg.astype(int)
        cls = {
            "active":                (~ig) & (ang == 0),
            "radial ghost":          ( ig) & (ang == 0),
            "angular face ghost":    (~ig) & (ang == 1),
            "angular corner ghost":  (~ig) & (ang == 2),
            "mixed radial x angular":( ig) & (ang >= 1),
        }
        rel = np.abs(da-db)/np.maximum(np.abs(da), 1e-300)
        relv = rel[neq]
        print(f"  {v}: {nd} cells differ (gz widths i,j,k = {gi},{gj},{gk}), "
              f"max rel {np.nanmax(relv):.2e}, median rel {np.median(relv):.1e}, "
              f"i range ({idx[:,3].min()}, {idx[:,3].max()})")
        for cname, sel in cls.items():
            n = int(sel.sum())
            if n == 0: continue
            print(f"      {cname:24s}: {n:8d}  max rel {np.nanmax(relv[sel]):.2e}  "
                  f"i range ({idx[sel,3].min()}, {idx[sel,3].max()})  "
                  f"j range ({idx[sel,2].min()}, {idx[sel,2].max()})  "
                  f"k range ({idx[sel,1].min()}, {idx[sel,1].max()})")
    return same

arms = [d for d in sorted(os.listdir('.')) if os.path.isdir(os.path.join(d,'E'))]
verdict = {}
for a in arms:
    print(f"--- ARM {a}")
    fp, P  = last(a+'/P')
    fe, E  = last(a+'/E')
    f0, F0 = last(a+'/F0')
    ff, F  = last(a+'/F')
    v = {}
    if P is None or F0 is None:
        print(f"  55800: P dump {fp}, F0 dump {f0}"); v['55800'] = "NO DUMP"
    else:
        print(f"  [55800]  P {fp} t={P['time']}  |  F0 {f0} t={F0['time']}")
        if P['time'] != F0['time']: print("  !! different times at 55800")
        v['55800'] = "IDENTICAL" if report('55800', P, F0) else "DIFFER"
    if E is None or F is None:
        print(f"  55801: E dump {fe}, F dump {ff}"); v['55801'] = "NO DUMP"
    else:
        print(f"  [55801]  E {fe} t={E['time']}  |  F {ff} t={F['time']}")
        if E['time'] != F['time']: print("  !! different times at 55801")
        v['55801'] = "IDENTICAL" if report('55801', E, F) else "DIFFER"
    verdict[a] = v
    print(f"ARM {a}: 55800 {v['55800']}   55801 {v['55801']}")
print("=============== SUMMARY")
print(f"{'arm':16s} {'P vs F0 @55800':>16s} {'E vs F @55801':>16s}")
for a in arms:
    print(f"{a:16s} {verdict[a]['55800']:>16s} {verdict[a]['55801']:>16s}")
print("""
READING RULE
  * 55800 DIFFER, and the differing cells are GHOSTS only  -> the restart's ghost fill is
    not the running ghost fill.  The class printed (radial / angular face / angular corner
    / mixed) names WHICH fill, and that is the answer.
  * 55800 IDENTICAL everywhere including ghosts -> the init path reproduces the whole
    array, the E/F method is sound, and the divergence at 55801 is made DURING the cycle:
    look at the arms and at the dt lines above.
  * any arm whose 55801 turns IDENTICAL names the state the restart does not restore.
  * onerank: 55801 IDENTICAL with 1 rank but DIFFER with 2 -> it is the MPI seam /
    coarse-buffer exchange state, not the physics.""")
PY
