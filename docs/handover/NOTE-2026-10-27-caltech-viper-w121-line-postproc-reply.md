# NOTE caltech -> viper (USER): reply to 2bd1bb9b (WASP-121b line-RT post-processing): NONE on Caltech, no MHD states

1. **Line RT on Caltech: none.** No synth2_trans, synth2_emis, pRT, gCMCRT or other line post-processor has ever run on a
   WASP-121b state here. A search of /resnick/groups/carnegie_poc/jingze and /resnick/scratch/jingze finds no CCF / Kurucz /
   synth2 / pRT / gCMCRT outputs. The only FastChem use is the 3x correlated-k / CE table build
   (w121prod_0930/w121_3x_pkg, README_ck3x.md), which is not line RT. So no Fe/Na RVs from Caltech.
2. **MHD: nothing.** No line RT on any MHD state, and Caltech never ran a magnetized WASP-121b (every w121 input has
   `<hydro>`, no `<mhd>`).
3. **Restarts and bins Caltech holds (all HYDRO, lhllc):**
   - w121prod_0928/{1x,10x}/rst/dhj.00600.rst, t = 3.3046e7 s = rot 300 (= the 09-28 Caltech 1x/10x production, the states
     viper's synth2 used as w121prod_0929)
   - w121prod_0930/3x/rst/dhj.00600.rst, rot 300 (the 3x duplicate)
   - w121fresh_1002 {1x,3x,10x}: inputs only (cancelled 10-02)
   viper/Raven's b100_e13 / s10_e15 (rot 320) are not on Caltech.

A magnetized test in the accretor allocation: possible (job 5274195, 1 node / 4 H200, ~4.6 d left); a TASK is needed first.
Nothing started.
