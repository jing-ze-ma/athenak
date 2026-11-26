# NOTE 2026-11-27zj caltech -> orion: restart-key warning check on 4 stars (your 27zj, RECONFATAL)

## Build and gate (worktree wt_orion1121/rcf at c9805a084, CUDA 12.9 + g++13, star_vet, node 02-09)
- BUILD rcf rc=0 wall=347 s TUs=280
- BUILD_OK rcf md5=ff182b014158b85243fb9019668ef9ef prov=c9805a08414d9edee17a4397778b18c3f5f95852
- GATE PASS rcf: athena_rcf_c9805a084 md5 ff182b01 node hpc-sm-02-09 pass1=3.720e-11 log date 11-25 22:55:26 end: Terminating on cycle limit
- prov: c9805a08414d9edee17a4397778b18c3f5f95852

## rstkeys.py per star (02-09, 4 GPUs, mpirun -np 4 inside a 4-GPU step; READ-ONLY hard links of the latest production rsts; decks = copies in reconfatal_1127/ of the run's own deck with its command-line keys put in; extra keys: the output blocks off + time/tlim=1e9 because every rst is at its run's tlim; the rsts the runs wrote are deleted)
### m35mms
- rst:  
```
srun: error: Unable to create step for job 7333072: Memory required by task is not available
```
### m35zams
- rst:  
```
srun: error: Unable to create step for job 7333072: Memory required by task is not available
```
### t9l62
- rst:  
```
srun: error: Unable to create step for job 7333072: Memory required by task is not available
```
### t9mlt
- rst:  
```
srun: error: Unable to create step for job 7333072: Memory required by task is not available
```

## rstkeys.py per star (02-09, 4 GPUs, mpirun -np 4 inside a 4-GPU step; READ-ONLY hard links of the latest production rsts; decks = copies in reconfatal_1127/ of the run's own deck with its command-line keys put in; extra keys: the output blocks off + time/tlim=1e9 because every rst is at its run's tlim; the rsts the runs wrote are deleted)
### m35mms
- rst: 898e74eb565f99b7018a37132c018b21 /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/m35mms.rst
```
a_deck rc 131 CLEAN WARNING 0
b_named rc 131 CLEAN WARNING 0
c_ronly rc 131 CLEAN WARNING 0
ROW deck_m35mms.athinput | warnings with DECK 0 | with DECK.named 0 | b_named vs c_ronly: 0 files, maxdiff 0
```
### m35zams
- rst: 863a98e8c3503466212e8df06cda4951 /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/m35zams.rst
```
a_deck rc 131 CLEAN WARNING 0
b_named rc 131 CLEAN WARNING 0
c_ronly rc 131 CLEAN WARNING 0
ROW deck_m35zams.athinput | warnings with DECK 0 | with DECK.named 0 | b_named vs c_ronly: 0 files, maxdiff 0
```
### t9l62
- rst: 66766f0240ae0dc1ad1587ebf10d01aa /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/t9l62.rst
```
a_deck rc 131 CLEAN WARNING 0
b_named rc 131 CLEAN WARNING 0
c_ronly rc 131 CLEAN WARNING 0
ROW deck_t9l62.athinput | warnings with DECK 0 | with DECK.named 0 | b_named vs c_ronly: 0 files, maxdiff 0
```
### t9mlt
- rst: 29a18958bea467d668dc6cb849c88a57 /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/t9mlt.rst
```
a_deck rc 131 CLEAN WARNING 0
b_named rc 131 CLEAN WARNING 0
c_ronly rc 131 CLEAN WARNING 0
ROW deck_t9mlt.athinput | warnings with DECK 0 | with DECK.named 0 | b_named vs c_ronly: 0 files, maxdiff 0
```

## VOID: the two "rstkeys rows" posts above (22fcb5584 1b3e763cd )
- Both posts are void. No run in them started the code.
  - First pass: srun refused the step (600G asked, the holder gives 480G per node).
  - Second pass: mpirun could not bind 4 ranks inside the step, so every run exited with rc 131.
  - Their "warnings 0 / 0 files" rows mean nothing.
- The third pass (mpirun --bind-to none) is running on 02-09. Its rows replace them.
- The build and tiny gate post stays valid: rcf c9805a084, md5 ff182b01, pass1 3.720e-11.

## rstkeys.py per star (02-09, 4 GPUs, mpirun -np 4 inside a 4-GPU step; READ-ONLY hard links of the latest production rsts; decks = copies in reconfatal_1127/ of the run's own deck with its command-line keys put in; extra keys: the output blocks off + time/tlim=1e9 because every rst is at its run's tlim; the rsts the runs wrote are deleted)
### m35mms
- rst: 898e74eb565f99b7018a37132c018b21 /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/m35mms.rst
```
a_deck rc 0 CLEAN WARNING 2
  stored  problem/mlt_src_time = start   (code default stage)
  stored  problem/star_sponge_form = lawson   (code default exp)
b_named rc 0 CLEAN WARNING 0
c_ronly rc 0 CLEAN WARNING 0
ROW deck_m35mms.athinput | warnings with DECK 2 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0
```
### m35zams
- rst: 863a98e8c3503466212e8df06cda4951 /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/m35zams.rst
```
a_deck rc 0 CLEAN WARNING 2
  stored  problem/mlt_src_time = start   (code default stage)
  stored  problem/star_sponge_form = lawson   (code default exp)
b_named rc 0 CLEAN WARNING 0
c_ronly rc 0 CLEAN WARNING 0
ROW deck_m35zams.athinput | warnings with DECK 2 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0
```
### t9l62
- rst: 66766f0240ae0dc1ad1587ebf10d01aa /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/t9l62.rst
```
a_deck rc 0 CLEAN WARNING 2
  stored  problem/mlt_src_time = start   (code default stage)
  stored  problem/star_sponge_form = lawson   (code default exp)
b_named rc 0 CLEAN WARNING 0
c_ronly rc 0 CLEAN WARNING 0
ROW deck_t9l62.athinput | warnings with DECK 2 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0
```
### t9mlt
- rst: 29a18958bea467d668dc6cb849c88a57 /resnick/groups/carnegie_poc/jingze/reconfatal_1127/rst_in/t9mlt.rst
```
a_deck rc 0 CLEAN WARNING 1
  stored  problem/star_sponge_form = lawson   (code default exp)
b_named rc 0 CLEAN WARNING 0
c_ronly rc 0 CLEAN WARNING 0
ROW deck_t9mlt.athinput | warnings with DECK 1 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0
```

## Valid rows (third pass; replaces the two void posts)
- The binary is rcf c9805a084 (md5 ff182b01, tiny gate 3.720e-11).
- Run on 02-09: 4 GPUs, mpirun -np 4 --bind-to none inside a 4-GPU step, READ-ONLY hard links of the latest production rsts.
- Decks are copies of each run's own production deck with its command-line keys put in (reconfatal_1127/deck_<star>.athinput).
- Extra command-line keys: output blocks off, and time/tlim=1e9 because every rst is at its run's tlim. rsts the runs wrote were deleted.
### m35mms
- rst m35mms.rst (hard link of the production file), md5 898e74eb565f99b7018a37132c018b21
- `ROW deck_m35mms.athinput | warnings with DECK 2 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0`
- `stored  problem/mlt_src_time = start   (code default stage)`
- `stored  problem/star_sponge_form = lawson   (code default exp)`
### m35zams
- rst m35zams.rst (hard link of the production file), md5 863a98e8c3503466212e8df06cda4951
- `ROW deck_m35zams.athinput | warnings with DECK 2 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0`
- `stored  problem/mlt_src_time = start   (code default stage)`
- `stored  problem/star_sponge_form = lawson   (code default exp)`
### t9l62
- rst t9l62.rst (hard link of the production file), md5 66766f0240ae0dc1ad1587ebf10d01aa
- `ROW deck_t9l62.athinput | warnings with DECK 2 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0`
- `stored  problem/mlt_src_time = start   (code default stage)`
- `stored  problem/star_sponge_form = lawson   (code default exp)`
### t9mlt
- rst t9mlt.rst (hard link of the production file), md5 29a18958bea467d668dc6cb849c88a57
- `ROW deck_t9mlt.athinput | warnings with DECK 1 | with DECK.named 0 | b_named vs c_ronly: 6 files, maxdiff 0`
- `stored  problem/star_sponge_form = lawson   (code default exp)`

## Verdicts
- `problem/star_sponge_form = lawson`: *intended*. It is in your list with the same value. It is stored in all four stars.
- `problem/mlt_src_time = start` (code default stage): *EXTRA*. It is stored in m35mms, m35zams and t9l62.
  - In those rsts the parameter block says `mlt_src_time = start   # Default value added at run time`. That makes it the start-up default of the code that wrote the rst (the production binary c4speed d979d5268), not a deck or hand-set key.
  - t9mlt names it in its own deck (`mlt_src_time = start   # T9MLT`), so t9mlt prints no warning for it. For t9mlt it is a hand-set key.
  - The user decides: name start (keeps the runs as they ran) or take the release default stage.
- No DIFFERS key: none of the 21 listed keys is stored at a value other than the list's.
  - The list's other 20 keys print nothing, because the stored value equals the code default of c9805a084.
- No reconfatal FATAL in any of the 12 runs. b_named equals c_ronly bitwise for every star (6 files, maxdiff 0).

## Item 3 (m35zams, t9mlt)
- Production decks, kept as DECK.named, which names exactly the listed keys at the stored values:
  - docs/handover/caltech_results/RECONFATAL/deck_m35zams.athinput.named: adds mlt_src_time = start and star_sponge_form = lawson.
  - docs/handover/caltech_results/RECONFATAL/deck_t9mlt.athinput.named: adds star_sponge_form = lawson.
- Both are copies; nothing under prod_* was changed. The rstkeys outputs of all 4 stars are in docs/handover/caltech_results/RECONFATAL/.

## Item 4 (m35mms, t9l62)
- The printed list is not equal to your list. It holds star_sponge_form (in the list) plus the EXTRA mlt_src_time = start.
- Of your list, only star_sponge_form differs from the c9805a084 code default. The other 20 are proven by silence: no warning means stored = code default = the template value.
- So the c9805a084 templates fit m35mms and t9l62 once mlt_src_time is decided.
