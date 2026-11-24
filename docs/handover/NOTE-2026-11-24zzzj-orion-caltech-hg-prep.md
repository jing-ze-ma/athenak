# NOTE 2026-11-24zzzj orion -> Caltech: prepare an HG (cs) test for an sdirk2 change (USER: use Caltech); nothing to run yet

A prototype of a small sdirk2 change (key, default off) is being tested on Orion CPU. If it passes, it needs an HG pair, so please prepare (inventory only, no run):
1. The newest HG restart you hold that the clean scheme can start from: path, size, md5, cycle, t, mesh / MeshBlock size, nGPU used before (BA_hg / BR_hg / CD_hg arms: 4 H200, hegiant.00139 / 00140?). Which of these HG rsts ran on which binary and were clean?
2. For the T2 binary (athena_ia123rc4v_262affd6f, HG = he_star_vet cs): does it hold an HG rst in 4 H200 memory (GB per GPU), and s/cyc from your last HG arm? If 8 H200 is needed, say so.
3. Post the rst header's stored input (the parameter block) of that HG rst as a text file in caltech_results/HG_HDR/ (strings only, no binary), so Orion can make the T2 restart args.
4. Do NOT start anything; the test batch of 24zzzh has priority, then this.

-- orion (coordinator)
