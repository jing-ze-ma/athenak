#!/bin/bash
# ridort-1014 (Raven worker RID): FULL IDORT for option R (arm r) and ROL2U (arm o, R + implicit_opac_lag=pred passes 2 + gas_newton_u4)
# for option R = rad_m1/implicit_tol_rel=0.1 + implicit_tol_rel_max=1e-4 (relative Picard stop)
# on fin51-r-fix b16b66db (straight runs bitwise to fin51-r 1bd0c2ce; default-off == FIN5.1 21c4839d). ONE binary, two arms per dir, concurrently, same threads:
#   ref = deck AS IS,  r = deck + the two keys (setkeys.py).  Outputs <OUT>/runs/<dir>/{ref,r} + <dir>/cmp.txt (one.sh).
# bash run_ridort.sh <plain CPU build of fin51-r> <cores of the node> <list> [<list> ...]
#   lists: list_a.txt (214 part-A = ENTHBE dirs), list_b.txt (35 = the families of the 20-dir CPU set). "<dir> <threads> <est s>".
#   Each dir takes 2 x threads cores; longest first. Edit the threads column (e.g. 32 for cloud / rs n1024 as in ENTHBE) if you like.
# env (defaults = Caltech): DECKS = the CVALL/wavefix tree with <dir>/in.athinput + args.txt [+ two.txt]; SHOCK = its shockref;
#   OUT = results root; PY = a python with numpy; BINCONV = <athenak>/vis/python (bin_convert.py, for the bin maxrel / reval).
HERE=$(cd $(dirname $0) && pwd); export BIN=$(readlink -f ${1:?}) HERE; C=${2:?}; shift 2
export DECKS=${DECKS:-/resnick/groups/carnegie_poc/jingze/wavefix_1013/caltech_wavefix/runs}
export SHOCK=${SHOCK:-/resnick/groups/carnegie_poc/jingze/wavefix_1013/caltech_wavefix/shockref}
export OUT=${OUT:-/resnick/groups/carnegie_poc/jingze/ridort_1014}; export PY=${PY:-python3} BINCONV=${BINCONV:-}
export ARMS=${ARMS:-ref r o} DECKS2=${DECKS2:-$HERE/decks_rk}
mkdir -p $OUT; md5sum $BIN | tee $OUT/md5.txt
for L in "$@"; do $PY $HERE/sched.py $L $C $HERE/one.sh | tee -a $OUT/rid_cmp_$(basename $L .txt).txt; done
echo "rc!=0: $(cat $OUT/rid_cmp_*.txt | grep -v 'rc=0/0' | wc -l)   BITWISE: $(cat $OUT/rid_cmp_*.txt | grep -o ': BITWISE' | wc -l)   DIFF: $(cat $OUT/rid_cmp_*.txt | grep -o ': DIFF' | wc -l)"
