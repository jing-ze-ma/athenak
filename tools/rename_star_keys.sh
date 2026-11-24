#!/bin/bash
# rename_star_keys.sh FILE ...: the <problem> keys he_*, hs_* -> star_* in place, in decks, key
# lists and launch scripts (he_gm = ..., problem/he_gm=...).  Values, paths and file names
# are kept; prints the number of renamed keys per file.  (The code also reads the old names.)
K='base_heat|bc|drive|esrc|gm|gravity|grid_dump|ic|nfine|opac|planck_table|seed|sponge|'
K+='ad_dead_skip|tau_int|wall|wb_grav_window|wind'
for f in "$@"; do
  perl -pi -e 'BEGIN {$n = 0}
    $n += s{(?:(?<=problem/)|(?<![\w./-]))h[es]_(?=('"$K"')(_\w+)?(?![\w/-]|\.\w))}{star_}g;
    END {print STDERR "$ARGV: $n keys\n"}' "$f"
done
