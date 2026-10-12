#!/bin/bash
# Render missing movie frames (Fig 3 phi=0 plane, Fig 4 r~36 Rsun shell) for every dump of
# the repro run (repro_4n, read only over the /raven/ptmp NFS mount) and re-assemble
# fig3.mp4 / fig4.mp4 (+ fig3.gif / fig4.gif only if <= 8 MB).  Rerun = only new dumps
# are reduced/rendered; colour ranges stay fixed (ranges.json; delete it to re-choose from
# the newest dump, then also delete frames/ so old frames are re-rendered).
set -e
M=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source /etc/profile.d/modules.sh 2>/dev/null || true
module load ${FFMPEG_MODULE:-ffmpeg/7.1} 2>/dev/null || true   # any ffmpeg on PATH
nice -n 10 python3 $M/movie_frames.py --arm ${1:-repro} > $M/last_frames.log
FPS=${FPS:-5}
for f in f3 f4 fint; do
  nice -n 10 ffmpeg -loglevel error -y -framerate $FPS -pattern_type glob \
    -i "$M/frames/$f/${f}_*.png" \
    -vf "scale='min(1600,iw)':-2:flags=lanczos,pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" \
    -c:v libx264 -preset slow -crf 23 -pix_fmt yuv420p -movflags +faststart \
    $M/fig${f#f}.mp4
  nice -n 10 ffmpeg -loglevel error -y -framerate $FPS -pattern_type glob \
    -i "$M/frames/$f/${f}_*.png" \
    -vf "scale=800:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse=dither=bayer" \
    $M/frames/fig${f#f}.gif.tmp.gif
  if [ $(stat -c %s $M/frames/fig${f#f}.gif.tmp.gif) -le 8000000 ]; then
    mv -f $M/frames/fig${f#f}.gif.tmp.gif $M/fig${f#f}.gif
  fi
done
ls -l $M/fig*.mp4 $M/fig*.gif 2>/dev/null || true
