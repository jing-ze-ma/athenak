#!/bin/bash
# fe_step.sh (fehook-1102, USER): fe_hires Fe I blueshift update for the dhj_r1 page, launched DETACHED by the catch-up step
# (own lock logs/fe.lock, so the page and the next refresh never wait for it). fe_update.sh picks a free TEST GPU itself (lock
# owner fehgpu; never production GPUs, never 01-05/01-09) or does a 1-dump CPU fallback. Then: images -> img/, html snippet
# (fe_section.py) -> logs/fe_section.html, quick page rebuild in the allocation, READY rewritten with the image list.
P=/resnick/home/jingze/ATHENAK/pages/dhj_r1; F=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast
FPY=/resnick/groups/carnegie_poc/jingze/fewind_1102/venv_fehgpu/bin/python
JOB=$1; NODE=$2
echo "fe start $(date +%T)" > $P/logs/fe.out
mkdir -p $P/logs/fe_img
bash $F/fe_update.sh --max 8 --img $P/logs/fe_img >> $P/logs/fe.out 2>&1
for f in fe_vs_rot.png fehires_vs_seidel.png; do [ -s $P/logs/fe_img/$f ] && cp $P/logs/fe_img/$f $P/img/$f.tmp && mv $P/img/$f.tmp $P/img/$f; done
$FPY $F/fe_section.py --img-prefix img/ > $P/logs/fe_section.html.tmp 2>> $P/logs/fe.out && [ -s $P/logs/fe_section.html.tmp ] && mv $P/logs/fe_section.html.tmp $P/logs/fe_section.html
echo "fe images + section $(date +%T)" >> $P/logs/fe.out
QH='module load python/3.11.6-gcc-11.3.1-rphh4kv py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui >/dev/null 2>&1; cd '$P'/tools && timeout 120 python3 ana.py html'
timeout 300 srun --overlap --jobid=$JOB -w $NODE -N1 -n1 -c 2 --mem=16G --gres=none -J dhjfe nice -n 10 bash -c "$QH" >> $P/logs/fe.out 2>&1 < /dev/null
H=$(grep '^dhj_r1 READY' $P/logs/fe.out | tail -1 | sed 's/^dhj_r1 READY //'); [ -n "$H" ] || H=$(head -1 $P/READY 2>/dev/null | sed 's/^READY [0-9-]* [0-9:]* //')
{ echo "READY $(date +%m-%d\ %H:%M) $H"; cd $P; echo index.html; ls img/*.png; } > $P/READY.tmp && mv $P/READY.tmp $P/READY
echo "fe done $(date +%T)" >> $P/logs/fe.out
