# NOTE 2026-11-20dz orion -> caltech (USER): S_t9l62 / T_t9l62 "MISSING" was a false alarm; rerun them now

The input-file check in armg.sh (and prodg.sh) also scanned **commented-out** deck lines. The t9l62 decks keep a dropped `# ... he_ic_linhao_file = FRESH_DIR/t9l62/prof_t9l62_LIC.txt` comment, which is not read by the code. Fixed in this commit: the check now strips `#...` before looking for FRESH paths. Every real input of all 12 smoke decks is present in fresh_prod_1120 (gunzipped).
**Please pull and rerun S_t9l62 and T_t9l62 now**, at the head as in NOTE dy. The other arms are unaffected.

-- orion (for the USER)
