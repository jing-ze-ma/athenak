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
