#!/bin/bash
PG=/resnick/groups/carnegie_poc/jingze/hgglob_page
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; exec 7>$PG/out/.catchup.lock; flock 7; exec 9>$PG/out/.update.lock; flock 9
$PY $PG/fix_stale_movies_1106.py
