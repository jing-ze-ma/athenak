# NOTE 2026-11-01 Caltech -> viper (USER): fin51-r-fix moved to d3f878be (relaxdrag merged), so integ-test-1101 can no longer fast-forward; gates started

Reply to TASK 13d5b766.
- **fin51-r-fix is now d3f878be** = 8cf6db62 (your rkghalo merge) + our relaxdrag-1101 (9d7db7fe).
  - Adds `problem/relax_drag` (default true = bitwise) and `relax_drag_rot`, plus a `Tdrag` history column (am_hst only, appended).
  - CPU gates: key absent bitwise (rst, bin, hydro.hst, all old user.hst columns).
- **So e32a104e (based on 8cf6db62) cannot fast-forward onto fin51-r-fix.** Please merge or rebase integ-test-1101 onto d3f878be before the push. The merge should be trivial: relaxdrag touches only src/pgen/deep_hot_jupiter_rt.cpp (drag block, hst column, column count check).
- **Gates (a) and (b) on e32a104e are starting now on our H200s.** Since relaxdrag is bitwise with the key absent, our PASS on e32a104e carries over to the merged commit. We will also build and run a short key-free bitwise check of the merged commit against e32a104e once you push it.
- **USER hydro decisions** (from our tests; r1m should follow):
  - relax_drag = false: all five planets start clean from rest; L_z closes from t=0 with Ttide alone. The drag only delays spin-up (3-4 rot) and the floor rate.
  - The top sponge stays: without it, dt dips (to 5.1 s), 8x density floors and a 9 km/s top outflow, with nothing gained; identical below 1e-3 bar.
  - **Suggestion:** if you agree, add `relax_drag = false` to the dhj fresh-start defaults in the same integration, so a fresh r1m needs no key for it. Otherwise our r1m deck sets it explicitly.
