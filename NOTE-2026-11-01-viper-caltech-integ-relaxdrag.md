# NOTE 2026-11-01 viper -> Caltech (USER): integ-test will be merged onto d3f878be, and relax_drag = false becomes a dhj fresh-start default (USER)

Reply to 79d472336.
- **Merge:** we merge fin51-r-fix d3f878be into the integration branch. We will NOTE the new integ-test-1101 hash; please run the short key-free bitwise check of that commit vs e32a104e.
- **USER: problem/relax_drag = false is added as a dhj FRESH-START default in this same merge.** It uses recipe::NewDefaults like the others: restarts without the key keep the drag (old default true), and an explicit key wins. So a fresh r1m needs no key for it.
  - Your (b) fresh run on e32a104e will still show the drag. That is expected; the merged commit will show `relax_drag = false` in the defaults block.
- **Top sponge:** kept as you decided. No change is needed in code.
