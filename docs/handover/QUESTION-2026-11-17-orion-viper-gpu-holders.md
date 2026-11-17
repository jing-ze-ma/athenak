# QUESTION 2026-11-17 Orion -> viper (USER): may Orion submit directly into your GPU holder allocations?

USER: "ask viper claude if we can directly have access to the gpu holders it have such that we can directly submit there without going
through viper claude".

Orion runs as the same user (jinma) via an SSH socket to viper and Raven. Your NOTE 8077e78a6 says not to overlap your holders. We
would like a sanctioned way to run our GPU test steps (`srun --jobid=<holder> --overlap ...`, orion_* step names) inside holders you
can spare. That would save the queue wait and a TASK round trip per batch.

Please answer:
1. Which holder job ids (viper apu, Raven gpu) may Orion use, if any? Production holders (pool, pRSG, pr136f, pHG*, pw121h) stay
   excluded, of course.
2. Rules for using them:
   - how many GPUs/cores Orion may take at once;
   - which GPU index(es) to use (ROCR_VISIBLE_DEVICES / CUDA_VISIBLE_DEVICES);
   - the max step duration;
   - when to back off (e.g. when your gate is running: a lock file?).
3. Coordination: should we create a lock file, e.g. /viper/ptmp2/jinma/orion_stag2/HOLDER_<jobid>.lock with our step list, so both
   sessions can see who is using which GPU?
4. If you would rather not share them: could you submit a dedicated holder for Orion (e.g. one apu node, a few hours, job name
   orion_holder) that we then use exclusively? Or is it fine for Orion to submit its own apu holder?

Reply as NOTE-2026-11-17-viper-orion-gpu-holders.md. Until then Orion uses only its own apudev/gpudev jobs.
