# TASK 2026-10-24 caltech -> viper (USER): move the 3-D accretor production (run12) to Caltech, 2 nodes

USER (10-24): "ask viper to move accretor run here on 2 nodes". Caltech = 2 nodes x 4 H200 = 8 GPUs (8 ranks, CUDA, -ffp-contract=off
equivalent = Caltech _nofma builds). RSG E and HG global keep running on Caltech; this run gets its own 2-node allocation.

## Please
1. Hand-over bundle of run12 (acc3d_prod_D_B1v6_extrap2_viper) in one directory on viper, for the USER to copy (no ssh caltech->viper):
   - the newest restart (+ md5), the deck exactly as run (incl. all command-line keys of the current link), the run.cfg / link script,
     hst/user files needed to continue the history (or say which to start fresh), and any input tables the deck reads.
   - binary info: commit (958cac780 r71 build?), configure line (problem gen, flags), so Caltech builds the same commit for CUDA.
   - the exact rsync command the USER should run to put it at /resnick/groups/carnegie_poc/jingze/acc3d_run12_move/.
2. Meshblocks: does the run's meshblock count split over 8 ranks? If not, say which rank count / decomposition Caltech should use.
3. Reference numbers for a port check: last 3 hst rows + dt / Picard mean of the last link, so Caltech can confirm a bitwise-or-
   round-off continuation in a short smoke before production.
4. Keep run12 running on viper until Caltech posts a NOTE that its continuation STARTED from your bundle; then stop viper's chain
   (USER decides if it differs).
5. Post a NOTE viper->caltech with the bundle path, md5s, rsync command and the answers above.
