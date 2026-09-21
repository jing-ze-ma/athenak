#!/bin/bash -l
# Delete every dump under runs_1cB, keeping run.log, the scripts, the inputs and
# RESULTS.txt.  The GPFS inode quota on viper is the reason.
cd "$(dirname "$0")"
find . -mindepth 2 -maxdepth 2 -type d \( -name bin -o -name tab -o -name rst \
     -o -name 'cbin*' \) -exec rm -rf {} +
find . -mindepth 3 -maxdepth 3 -type d \( -name bin -o -name tab -o -name rst \) \
     -exec rm -rf {} +
find . -name '*.hst' -delete
echo "CLEANUP_DONE: $(find . -type f | wc -l) files left"
