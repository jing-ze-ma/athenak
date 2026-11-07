#!/bin/bash
# 11-01 USER "make updating pages fast and easy". One command for all live pages:
#   pages.sh        -> quick text refresh of every page (login, ~1 min) + detached in-alloc catch-ups for new dumps
#   pages.sh status -> one line per page ("READY <time> t = ..." means movies are done: publish index.html + files listed in <page>/READY)
J=/resnick/groups/carnegie_poc/jingze
for n in rsgE hggl; do bash $J/quickpage.sh $n $1; done
[ -x /resnick/home/jingze/ATHENAK/pages/dhj_r1/refresh.sh ] && { [ "$1" = status ] && echo "dhj_r1: $(head -1 /resnick/home/jingze/ATHENAK/pages/dhj_r1/READY 2>/dev/null || echo 'not ready')" || bash /resnick/home/jingze/ATHENAK/pages/dhj_r1/refresh.sh; }
