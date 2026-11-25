#!/bin/bash
# mov.sh STAR: fixed limits over the whole run, all frames, one H.264 movie (<= 8 MB) + poster
cd /resnick/groups/carnegie_poc/jingze/docfig_1125; unset PYTHONPATH PYTHONHOME; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; PY=/resnick/groups/carnegie_poc/jingze/palette_1125/venv_pal/bin/python; s=$1; FF=/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin
$PY docfig.py $s limits | tail -1 | cut -c1-200; NPROC=${NPROC:-10} $PY docfig.py $s frames | tail -1
$PY - <<P
import glob,os
D='frames/$s/'; F=sorted(glob.glob(D+'f_*.png')); bad=[f for f in F if open(f,'rb').read()[-12:].find(b'IEND')<0]; print('frames',len(F),'bad',len(bad)); n=len(F); st=max(1,round(n/280)); pick=F[::st]
os.makedirs('seq/$s',exist_ok=True)
for q in glob.glob('seq/$s/*.png'): os.remove(q)
for i,f in enumerate(pick): os.symlink(os.path.abspath(f),'seq/$s/%05d.png'%i)
import json; json.dump(dict(n=n,step=st,used=len(pick)),open('work/mov_$s.json','w')); print('used',len(pick),'of',n)
P
N=$(ls seq/$s | wc -l); FR=$($PY -c "n=$N; print(min(30,max(n/10.0,1)))")
for crf in 22 26 30; do $FF/ffmpeg -y -loglevel error -framerate $FR -i seq/$s/%05d.png -vf "scale=1200:-2:flags=lanczos,fps=24" -c:v libx264 -preset slow -crf $crf -pix_fmt yuv420p -movflags +faststart out/mov_$s.mp4; SZ=$(stat -c %s out/mov_$s.mp4); [ $SZ -le 4000000 ] && break; done
F0=$(ls frames/$s/*.png | tail -1); cp $F0 out/mov_${s}_poster.png; $FF/ffprobe -v error -count_frames -select_streams v -show_entries stream=nb_read_frames,width,height,r_frame_rate,duration -of csv=p=0 out/mov_$s.mp4; ls -la out/mov_$s.mp4 | awk '{print $5}'
