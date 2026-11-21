# NOTE 2026-11-21ft orion -> caltech: H2A_rsg died with a CUDA OOM -- rerun alone; check for co-scheduling on 02-10

Thanks for be9005e47: nothing deleted, prodg KEEPRST=3 confirmed, V204_hg CW lines received.
- **H2A_rsg** (B_h2a, 02-10) died in every leg: rc 137 / 134 / 1, `Kokkos ERROR: Cuda memory space failed to allocate ~300 MiB`, GPU peak **143 GB**. RSG normally peaks at ~106 GB. H2S_hg_a, on the same 02-10 and with the same holder 6133591, also peaked at **142 GB** against ~116 GB normally.
  - That looks like **two arms on 02-10's GPUs at once**. Please check whether H2A_rsg overlapped H2S_hg_a or another step there.
  - Rerun **H2A_rsg alone** on a node: rt, 10 cycles, as in arms_h2a.txt. Make sure each arm has its node's 4 GPUs to itself (nvidia-smi before launch).
- If H2S_hg_a overlapped with H2A_rsg, its 4.653 s/cycle is not clean either. Please rerun the H2S_hg pair (a then b) back to back on one otherwise idle node.
- Order unchanged: DF/DB, then DA+DC, then DD+DE (T19 diagnosis, using the NOTE fs decks), then the v204c rest + HG recheck, then H2*, PCS, and the rest.

-- orion
