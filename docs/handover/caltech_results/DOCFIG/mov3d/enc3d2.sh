#!/bin/bash
# enc3d.sh: 1080p movies of the five stars (existing frames)
cd /resnick/groups/carnegie_poc/jingze/mov3d_1125; mkdir -p out; FF=/central/software9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.3.1/ffmpeg-6.0-eim3de2umcr2h3ugaxw7as4pplxya727/bin
X="-c:v libx264 -preset slow -pix_fmt yuv420p -movflags +faststart"
plain(){ s=$1; fr=$2; out=$3; for crf in 20 22 24 26 28 30; do $FF/ffmpeg -y -loglevel error -framerate $fr -start_number ${START:-0} -i frames/${s}_%05d.png -vf fps=$out $X -crf $crf out/mov3d_$s.mp4; SZ=$(stat -c %s out/mov3d_$s.mp4); echo $s crf $crf $SZ; [ $SZ -le 6000000 ] && break; done; }
fade(){ s=$1; for crf in 22 24 26 28 30 32; do $FF/ffmpeg -y -loglevel error -framerate 12 -start_number 0 -i frames/${s}_%05d.png -filter_complex "[0:v]format=yuv420p,split[a][b];[a]minterpolate=fps=30:mi_mode=blend[m];[b]fps=30,crop=iw:118:0:0[t];[m][t]overlay=0:0:shortest=1" $X -crf $crf out/mov3d_$s.mp4; SZ=$(stat -c %s out/mov3d_$s.mp4); echo $s crf $crf $SZ; [ $SZ -le 8000000 ] && break; done; }
case $1 in bsg) plain bsg 12 24;; m35mms|m35zams) fade $1;; t9l62) START=76 plain t9l62 15 30;; t9mlt) plain t9mlt 15 30;; esac
cp $(ls frames/$1_*.png | tail -1) out/mov3d_$1_poster.png
$FF/ffprobe -v error -count_frames -select_streams v -show_entries stream=codec_name,pix_fmt,width,height,r_frame_rate,nb_read_frames:format=duration -of csv=p=0 out/mov3d_$1.mp4
