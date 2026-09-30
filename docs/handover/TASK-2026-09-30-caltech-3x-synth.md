# TASK 2026-09-30 (Caltech -> viper): what Caltech needs to run synth.py phase on the 3x run

**WITHDRAWN 09-30 (Caltech does it itself; no viper action needed).** origin/dhj-olr-dump merges cleanly.  Caltech built
793e03c3 + dhj-olr-dump (local commit 2f667e95, job 3668813), then olr dumps from the rot-254 restarts of all three arms
and the rot-300 restarts of 1x/10x (job 3668862), then synth.py phase.  Answers to questions 2-5 are still welcome as FYI.

Context: Caltech analysed 3x at rot 254 against 1x/10x with the pushed w121 pipeline
(docs/handover/scripts/w121, 74f78c6e): ana300, synth winds, deepmix.  Our 1x at rot 300 reproduces viper's
ana_rot300 numbers.  Results: /resnick/groups/carnegie_poc/jingze/w121prod_0930/ana_3x_rot254/RESULTS.md.
Short version: 3x has settled at depth, is sub-adiabatic, and has no 10x-like deep cooling; jet/KE/deep T lie
between 1x and 10x; dt 7.7 s, and 3.1x the 1x T-floor cells (364 at p > 1e-6 bar).
**`synth.py phase` has not been run** (phase offset, NIRSpec night flux, eps, T_day/T_night).  It needs the
problem/olr_dump output.  That key exists only on origin/dhj-olr-dump (cd870c6f..2efb2851), which is not merged
into rt-integration.

Questions for viper (please answer in this file or in a NOTE, then push):
1. **Recipe.** How exactly did you produce the olr files for synth_rot300 / synth_rot300_10x?
   - Which restarts: how many snapshots, and over which rotation window?
   - Run settings: the nlim you used after the restart, and how you made the last cycle contain a full ck call
     (ck_impl_every 4: "reads the arrays of the LAST ck call").
   - Any other keys changed.
   - CPU or GPU, and the rank count.
   Please paste the athinput diff and the command line.
2. **Branch.** Is dhj-olr-dump safe to merge into rt-integration (default off, bitwise)?  If it is, please merge it.
   Otherwise, is Caltech OK to build `rt-integration + dhj-olr-dump` locally for these dumps only?
3. **EOS table.** synth.py (and dhjcs / ana300) read /viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt.
   Which key or tool writes it, and is it metallicity-dependent?  Caltech's stand-in is eos_dump_dhj.cpp, built from
   eos_composition.hpp with each run's X/Y/[M/H]; it is in the results dir above.  Please also push dclib.py (we used
   a stand-in).
4. **3x specifics.** Does synth.py phase need anything metallicity-specific beyond RTOP/R0/R1/CSTR (for example band
   edges taken from the dump header, or the stellar spectrum)?
5. **Which 1x/10x runs?** Were your synth_rot300 phase results computed from the Caltech 1x/10x restarts, or from
   viper's own w1x/w10x?  If from viper's, Caltech will also dump 1x/10x at rot 254, so that all three arms are
   compared on the same footing.

Caltech plan once answered: 3x link 2 (3648795) reaches rot 300 on ~09-30 20:00 PDT.  Then take olr dumps from the
3x restarts (and 1x/10x if needed) as a short 2-H200 job, run synth.py phase, and add the output to RESULTS.md.
