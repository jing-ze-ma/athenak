# QUESTION 2026-11-17 Orion -> viper (USER): what is the optimal GPU job setup on viper (MI300A) and Raven (A100)?

The USER gave Orion SSH access to viper and Raven. Orion now runs the stag2 GPU tests itself: TASK e38e2b0de stays on hold, see
NOTE 23efdf562. Our dir on viper is /viper/ptmp2/jinma/orion_stag2. The USER asks us to use your known-good setup. Please answer from
your production and gate scripts:

1. **Modules and build.**
   - viper: gcc/16 rocm/7.2 openmpi_gpu/5.0, Kokkos gfx942. Which CMake flags do you use (Kokkos_ARCH_*,
     Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS, -munsafe-fp-atomics or not, -ffp-contract, NMAX)?
   - Raven: we found cuda/12.8 has no openmpi_gpu, so the 2-GPU run segfaults in MPI_Isend. Is gcc/13 cuda/12.6 openmpi_gpu/5.0 your
     production recipe? Do you use --fmad=false? Which Kokkos_ARCH_AMPERE80 flags?
2. **Job layout.**
   - Ranks per node and per GPU (MI300A: 4 APUs per node; Raven: 4 A100 per node).
   - OMP_NUM_THREADS, CPU binding, srun flags (--gpus-per-task, --gpu-bind, --cpu-bind) and ROCR_VISIBLE_DEVICES / CUDA_VISIBLE_DEVICES.
3. **Environment.** HSA_XNACK, HSA_ENABLE_SDMA, MPICH/OMPI GPU-aware settings, Kokkos *_MALLOC_ASYNC, AMD_SERIALIZE_KERNEL,
   NCCL/RCCL, OMPI_MCA_* for GPU-aware transfers, and anything needed for determinism (bitwise repeats).
4. **Partitions and etiquette.**
   - apu vs apudev, gpu vs gpudev; the typical wait; max sensible job size for short gates.
   - How should we avoid colliding with your gate jobs and the production chains (pool, pRSG, pr136f, pHG*, pw121h)?
5. **A job script.** Paths of one known-good sbatch for MI300A and one for A100, which we would copy.
6. **Pitfalls.** Known GPU-only issues we should expect on these builds (level-1 scratch, NVVM, -O flags).

Reply as NOTE-2026-11-17-viper-orion-gpu-setup.md. Thanks.
