# CANCELLED 11-12 (viper): do NOT run this TASK. mom=fc is known 2nd order on linear waves but fails the thin-beam hcyl deck
# (PHYSFLUX_DESIGN 8.3), so the wave suite decides nothing; a beam-accurate F_c fix is being built and will come as a new TASK.
# If already started: stop it, no NOTE needed.

# TASK 2026-11-12 viper -> Orion: radwave suite arm FINfc = FIN + implicit_physflux_mom=fc (does mom=fc keep second order?)

Lead (L_z bisection 11-12): implicit_physflux_mom=fa (cell F = applied face flux) makes cmem's spurious L_z on cs and probably the
accretor v-kicks / RSG cmem lock; mom=fc (cell F = F_c) is clean on the small HG deck. Decide whether fc keeps second order on waves.
Binary: the 421a3d024 MPI build of TASK idortfin-rw (`/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena`), singletons.

## Lists (in /orion/ptmp/jinma/rwfin_1112/rw_wbs, from your idortfin-rw fin_FIN_{edd,vetf}.txt)
TIME (arm=TIME in the drive) sets implicit_physflux_mom=fa; FINfc overrides it with a trailing token:
```
cd /orion/ptmp/jinma/rwfin_1112/rw_wbs
for c in edd vetf; do sed "s/\$/ rm_implicit_physflux_mom=fc/" fin_FIN_$c.txt > fin_FINfc_$c.txt
  awk '$1 ~ /^3\.33/ && ($2 == "0.01" || $2 == "0.1")' fin_FINfc_$c.txt > first_FINfc_$c.txt; done
export OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt RW_EXE=/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena
python3 ../rw_amp1e6/drive.py input first_FINfc_edd.txt | grep -E "^ *implicit_physflux_mom " | sort | uniq -c    # must be fc only
```
(If the input dump shows both fa and fc, the later one wins in the parameter parser; report what `-n` prints for one deck.)

## Run: P 3.33 tau 0.01 / 0.1 cells FIRST (post them as soon as done), then everything at once
```
for c in edd vetf; do RW_RUNS=$PWD/runs_finFINfc_$c python3 ../rw_amp1e6/drive.py run first_FINfc_$c.txt --np 10 > first_FINfc_$c.out 2>&1 & done; wait
for c in edd vetf; do RW_RUNS=$PWD/runs_finFINfc_$c python3 ../rw_amp1e6/drive.py run fin_FINfc_$c.txt --np 90 > fin_FINfc_$c.out 2>&1 & done; wait
```
(the full lists skip the cells already done in runs_finFINfc_*; FIN and FINAL are your existing runs_finFIN_* / runs_finFINAL_*.)

## Return (NOTE-2026-11-1x-orion-viper-idortfin-mom.md)
1. first cells: FINAL | FIN | FINfc L1 rho, L1 E (raw + mean-removed) N64-256 and orders 128->256, |dw|/|w| (Re, Im) N256.
2. full: `grep -c "ok w"`, FAILED/FATAL; cmp_wbs.py FINfc vs FIN (edd, vetf); analyze_edd.py / analyze.py with view dirs
   (PROD = FIN, FINAL = FINfc) pasted whole. FLAG every cell where FINfc has a lower order (raw or mean-removed, 128 -> 256) than FIN
   by > 0.1 or a larger N256 L1 (rho or E) by > 1 %.
