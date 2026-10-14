#!/bin/bash
# mem-1014 (Raven worker MEM): fin4-mem (memory reductions, must be BITWISE) vs FIN4 on IDORT dirs, plain CPU, decks AS IS + args.txt.
# bash run_mem.sh <MEM = plain CPU fin4-mem build> <FIN4 = your plain CPU fin4-1014 d505d21f build> [parallel=4] [list=list_mem.txt]
#  per dir <MB>/runs/<dir>/mem/: O = FIN4, N = MEM, NB=10 cycles each (same threads, the list's count), final rst byte-equal
#  after <par_end> (cmprst.py) ; dirs in RS also: MEM NB/2 -> -r -> NB, final rst vs N (restart bitwise).
#  -> <dir>/mem/cmp.txt ; summary lines "MEM <BITWISE|DIFF|FAIL> <dir> [restart: ...]" on stdout.
HERE=$(cd $(dirname $0) && pwd); BN=$(readlink -f ${1:?}); BO=$(readlink -f ${2:?}); P=${3:-4}; LIST=${4:-$HERE/list_mem.txt}; NB=${NB:-10}
MB=${MB:-/resnick/groups/carnegie_poc/jingze/mbband_1013/caltech_mbband}; EB=${EB:-/resnick/groups/carnegie_poc/jingze/enthbe_1013/caltech_enthbe}
RS=${RS:-"sp/rsw_T/HFfZ_1841da_1 cart/cb/FIN3_r4 sp/core_L16_u/HFfZ_1841da_128"}
shk() { sed -i -e "s#/raven/ptmp/jinma/suite_1011/shock/#$1/shockref/#" -e "s#/viper/ptmp2/jinma/fin3_1012/shockref/#$1/shockref/#" $2; }
one() { r=$1; nt=${2:-4}; d=$MB/runs/$r; [ -f $d/in.athinput ] || d=$EB/runs/$r; [ -f $d/in.athinput ] || { echo "MEM NODIR $r"; return; }; m=$d/mem
  if [ ! -f $m/cmp.txt ]; then rm -rf $m; mkdir -p $m; cp $d/in.athinput $m/in.athinput; shk $MB $m/in.athinput
    grep -qE '^ *file_type *= *rst' $m/in.athinput || printf '\n<output99>\nfile_type = rst\ndt = 1.0e30\n' >> $m/in.athinput
    for x in O N; do b=$BN; [ $x = O ] && b=$BO
      ( mkdir -p $m/$x; cd $m/$x; OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $b -i ../in.athinput -d $m/$x $(cat $d/args.txt) time/nlim=$NB > run.log 2>&1; echo $? > rc.txt ) & done; wait
    fo=$(ls $m/O/rst/*.rst 2>/dev/null | sort | tail -1); fn=$(ls $m/N/rst/*.rst 2>/dev/null | sort | tail -1)
    if [ -n "$fo" ] && [ -n "$fn" ]; then s=$(python3 $HERE/cmprst.py $fo $fn | awk '{print $1}'); else s="FAIL(rc=$(cat $m/O/rc.txt)/$(cat $m/N/rc.txt))"; fi
    rs=""
    if echo " $RS " | grep -q " $r "; then
      ( mkdir -p $m/H; cd $m/H; OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -i ../in.athinput -d $m/H $(cat $d/args.txt) time/nlim=$((NB/2)) > run.log 2>&1 )
      ( mkdir -p $m/R; cd $m/R; OMP_NUM_THREADS=$nt OMP_PROC_BIND=false $BN -r $(ls $m/H/rst/*.rst | sort | tail -1) -d $m/R time/nlim=$NB > run.log 2>&1 )
      rs="restart: $(python3 $HERE/cmprst.py $fn $(ls $m/R/rst/*.rst 2>/dev/null | sort | tail -1) 2>&1 | awk '{print $1}')"
    fi
    echo "$s $r $rs" > $m/cmp.txt; rm -rf $m/O/rst $m/N/rst $m/H/rst $m/R/rst
  fi
  echo "MEM $(cat $m/cmp.txt)"; }
export -f one shk; export HERE BN BO NB MB EB RS
md5sum $BN $BO
grep -v "^#" $LIST | xargs -P $P -L 1 bash -c 'one "$@"' _ | tee $HERE/mem_cmp.txt
echo "COUNTS: $(awk '{print $2}' $HERE/mem_cmp.txt | sort | uniq -c | tr '\n' ' ')"
