#!/bin/bash
# pub.sh: stage the RSG E page for publishing with ONE cheap output. Rebuilds first unless NOBUILD=1 (autoupd.sh keeps it fresh).
# Prints: what changed in the page TEXT since the last staging (base64 stripped) + one READY line. Then publish (Claude):
#   Artifact publish file_path=/u/jma20/ATHENAK/pages/rsgE/rsgE_deltaai.html url=https://claude.ai/artifact/JegJZ1mkx7nqHqz3pk2aU9
#            files={"movies/rsgE_eq.mp4": ".../pages/rsgE/movies/rsgE_eq.mp4", "movies/rsgE_map090.mp4": ".../movies/rsgE_map090.mp4"}
H=/work/nvme/bivj/jma20/rsgE_page; S=/u/jma20/ATHENAK/pages/rsgE; L=$H/pub.log
if [ -z "$NOBUILD" ]; then flock $H/.lock nice -n 10 $H/update.sh > $L 2>&1 || { echo "BUILD FAIL: tail -20 $L"; tail -5 $L; exit 1; }; fi
P=$H/out/publish; [ -f $P/index.html ] || { echo "no page built"; exit 1; }
strip(){ sed -E 's/base64,[A-Za-z0-9+\/=]{100,}/base64,IMG/g' "$1" | sed -E 's/<[^>]+>/ /g; s/[[:space:]]+/ /g' | tr '.' '\n' | sed '/^ *$/d'; }
flock $H/.lock cp $P/index.html $S/rsgE_deltaai.html; cp $P/movies/*.mp4 $S/movies/
strip $S/rsgE_deltaai.html > $S/.text.new
if [ -f $S/.text.pub ]; then diff $S/.text.pub $S/.text.new | grep '^>' | cut -c1-200 | head -${MAXL:-25}; else echo "(first staging: full text)"; head -60 $S/.text.new; fi
mv $S/.text.new $S/.text.pub
echo "READY $(grep -o 'built [^(]*([^)]*) by [^<]*' $S/rsgE_deltaai.html | head -1 | cut -c1-120) | html $(du -k $S/rsgE_deltaai.html | cut -f1) kB, movies $(du -ck $S/movies/*.mp4 | tail -1 | cut -f1) kB"
