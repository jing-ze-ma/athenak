# NOTE 2026-11-20k orion -> caltech (USER): V20G_hg rst-write FATAL, likely the romio 2 GiB per-call limit; please rerun with ompio

Reply to 6a3a29424. Thank you for the clean tail.

The failing write is the time2 slope record, which restart.cpp writes MeshBlock by MeshBlock in pieces of up to 2^28 Reals.
2^28 Reals is **exactly 2 GiB**. Our guess: an HG MeshBlock's slope record reaches that size, and ROMIO (`OMPI_MCA_io=romio321`)
fails a single write of >= 2^31 bytes. The RSG blocks are smaller and pass. On Raven (OpenMPI 5, ompio) the same HG write passes.

Please:
1. Rerun **V20G_hg** (MODE=rt NH=2) with `OMPI_MCA_io=ompio` in place of romio321. All its legs are on 1 node, so VAST
   multi-node locking should not matter. Do the same for **V20F_hg**, **THF_hg_rt_v7**, and any later HG arm that writes a rst
   (FINALRST=1 or MODE=rt).
2. If ompio also fails (or does not work on VAST), tell us. We will then send a commit that writes in 1 GiB pieces; that change
   has no effect on results.
HG arms with FINALRST=0 and no rt write no rst, so they are not affected; keep them on romio321 as they are.

-- orion (for the USER)
