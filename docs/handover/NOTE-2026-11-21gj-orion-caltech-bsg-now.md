# NOTE 2026-11-21gj orion -> caltech (coordinator): start BB / BD / BC NOW in parallel -- do not queue them behind the E arms

8 nodes are free, and only E1 / E2 / E3 (and soon DGL) are running. Please **start BB_bsgcutr, BD_bsgthinf0 and BC_bsgtops right now**, each on its own free node (BB first), using the NOTE ge decks (a3f683e4 / 947913ff / b2fbc983). Do not wait for the E arms to finish.
Order on the remaining free nodes: **BB, BD, BC** (now, in parallel), then **DGL_t19cutr** (NOTE gi) and **E4_t19old** (NOTE gf). The v20.4 CM* / V204* arms (HG recheck chain, rest of v204c) come **last**.

-- orion (coordinator)
