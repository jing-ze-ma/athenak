# NOTE 2026-11-20l orion -> caltech (USER): v20g2 = fix for the HG rst write; 2 arms, top priority

New build **B_v20g2** = caltech-runs/v20g2, commit 240d9eb3258e7a156a6a6db6751926cb975165cb. It is v20g plus one commit: restart writes go out in 1 GiB pieces instead
of 2 GiB, and the read-ahead cap is 1 GiB per rank. The file bytes are unchanged. Same cmake line as before.
Arms: `docs/handover/caltech_runs_1121/arms_v20g2.txt`. Run them **with your DEFAULT io (OMPI_MCA_io=romio321)**, which is the
case that failed:
- V20G2_hg: MODE=rt NH=2, as V20G_hg.
- V20G2_rsg: MODE=rt NH=2. Its hst must equal V20G_rsg bitwise; please put that xcheck in the result line.
The ompio rerun of V20G_hg (NOTE 03aa43e0b) stays as a cross-check, after these two.

-- orion (for the USER)
