#!/bin/bash
# mov5.sh STAR: 2-D docs movie from frames23 (v_r at tau_R = 2/3): m35 12 dumps per s, 30 fps cross-faded; others as mov.sh
cd /resnick/groups/carnegie_poc/jingze/docfig_1125; s=$1; FF=/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin; mkdir -p out5; X="-c:v libx264 -preset slow -pix_fmt yuv420p -movflags +faststart"; N=$(ls frames23/$s/*.png | wc -l); FR=$(awk "BEGIN{x=$N/10.0; if(x>30)x=30; print x}")
case $s in
 m35*) for crf in 22 24 26 28 29 30 32; do $FF/ffmpeg -y -loglevel error -framerate 12 -start_number 0 -i frames23/$s/f_%05d.png -filter_complex "[0:v]scale=1200:-2:flags=lanczos,format=yuv420p,split[a][b];[a]minterpolate=fps=30:mi_mode=blend[m];[b]fps=30,crop=iw:42:0:0[t];[m][t]overlay=0:0:shortest=1" $X -crf $crf out5/mov_$s.mp4; SZ=$(stat -c %s out5/mov_$s.mp4); echo $s crf $crf $SZ; [ $SZ -le 6000000 ] && break; done;;
 *) for crf in 22 26 30; do $FF/ffmpeg -y -loglevel error -framerate $FR -start_number 0 -i frames23/$s/f_%05d.png -vf "scale=1200:-2:flags=lanczos,fps=24" $X -crf $crf out5/mov_$s.mp4; SZ=$(stat -c %s out5/mov_$s.mp4); echo $s N $N crf $crf $SZ; [ $SZ -le 4000000 ] && break; done;;
esac
cp $(ls frames23/$s/*.png | tail -1) out5/mov_${s}_poster.png; $FF/ffprobe -v error -count_frames -select_streams v -show_entries stream=nb_read_frames,width,height,r_frame_rate:format=duration -of csv=p=0 out5/mov_$s.mp4
