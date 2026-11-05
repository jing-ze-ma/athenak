# NOTE 2026-11-05b viper -> Caltech: GPU gap follow-up - Raven checked; fin51-r-fix fully in prod-1102

This is an addendum to NOTE-2026-11-05-viper-caltech-gpu-gap.md.

## Raven (socket re-opened by the user; checked 11-05)
- gpu partition:
  - nodes: 4 idle, 96 mixed, 84 allocated;
  - 508 jobs pending: 182 on priority, the rest on dependencies or held;
  - limits: MaxTime 24 h, MaxNodes 128.
- Our share: **mpa_gpu fair-share factor 0.44 with zero recent usage**, i.e. good priority, much better than Viper (0.08).
- /raven/ptmp is 97 % full overall. Our quota could not be read (no lfs/quota on the login node). Check the free space before copying ~160 GB.
- => **Raven is the machine to use for the gap.** The memory caveat stands: A100 40 GB per GPU, so please send the per-block GPU memory high-water mark for RSG E and HG.

## Branches
- `git cherry fork/prod-1102 fork/fin51-r-fix` is empty: every fin51-r-fix commit (head c38bd1981) is in prod-1102.
- prod-1102 is the production branch (USER). viper will commit no new production work to fin51-r-fix.

## Next (needs the user)
1. Which runs move.
2. The restart copy route (user-driven scp/rsync).
3. Then viper/Raven builds prod-1102 (8cd7f0a9f or newer) with build_inc_raven.sh and smokes the exact binary + deck before any chain.
