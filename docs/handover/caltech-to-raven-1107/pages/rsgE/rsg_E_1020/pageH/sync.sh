#!/bin/bash
# pageH/sync.sh (10-29): the RSG E page's view of the re-run: prodG dumps/hst BEFORE prodH's restart (1675 d) + all of prodH
# (prodH reuses prodG's dump numbers from 1676 d; prodG's own >=1676 d history is archived as *_prodG_post1675 in the page caches)
# pagefast-1102: links in one python process (only new ones), the prodG hst prefix cached (keyed by prodG hst size/mtime + t0);
# output byte-identical to sync.sh.pre_fast_1102 (was ~230 s, mostly per-file ln/sed forks).
R=/resnick/groups/carnegie_poc/jingze/rsg_E_1020; P=$R/pageH
/usr/bin/python3 $P/sync_links.py $R > /dev/null
for h in rsg.user.hst rsg.hydro.hst; do
  t0=$(grep -v -m1 '^#' $R/prodH/$h | awk '{print $1}')
  key="$t0 $(stat -c '%s %Y' $R/prodG/$h)"; C=$P/.prefix_$h
  if [ "$(cat $C.key 2>/dev/null)" != "$key" ]; then
    awk -v t0=$t0 '/^#/ || $1+0 < t0+0' $R/prodG/$h > $C.tmp && mv $C.tmp $C && echo "$key" > $C.key
  fi
  { cat $C; grep -v '^#' $R/prodH/$h; } > $P/$h.tmp && mv $P/$h.tmp $P/$h
done
ln -sfn $R/prodH/rsg.log $P/rsg.log
ln -sfn $R/prodH/rsg.x1grid.txt $P/rsg.x1grid.txt
