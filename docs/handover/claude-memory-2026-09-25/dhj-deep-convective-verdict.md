---
name: dhj-deep-convective-verdict
description: 09-25 deepconv_0925 -- dhj deep is convective below the RCB ~3 bar (nabla_rad/nabla_ad 7 at 10 bar, ~790 at 200); the 0.9 nabla_ad IC is mixed to neutral within ~60 rot and the deep COOLS 74-163 K to the cooled 3-bar entropy; fix = exact-adiabat IC on the relaxed adiabat
metadata:
  type: project
---
/viper/ptmp2/jinma/deepconv_0925/README.md (production EOS + conduction Rosseland table + pgen IC arrays dumped by an
analysis pgen). RCB 2.95-3.25 bar everywhere (day/night/terminator, hyd4 rot 102, prod4 rot 84). Radiation carries
only 0.1-15 % of F_int below it. IC: adjust_ad_pT_arr switches at 3.15 bar/3911 K to 0.9 nabla_ad. Both runs: deep
became nabla/nabla_ad 0.99 (3-100 bar; 0.95 at the wall), cooled 163/146/110/74 K at 10/30/100/200 bar (s -0.48..-0.09
k_B/m_u); a cooling front from ~3 bar reached 10 bar at rot ~3, 200 bar at rot ~20; settled by rot 60-68 (90 % by
46-56); last 20 rot <= 0.05 K/rot; hydro = MHD to 2 K. t_th = p c_p T/(gF) 4.6e4..6.1e5 rot at 10..200 bar -> the
drift is dynamical mixing, not radiative. The substellar picket-fence IC was too hot aloft (1 bar 3740 -> 3350 K).

**Why:** decides the dhj initial condition and whether MLT is needed.
**How to apply:** IC = exact adiabat below ~3 bar on the RELAXED adiabat (~3716 K at 3 bar), which the ck RCE global-mean
profile ([[sparc-sponge-campaign-0925]]) should give; MLT not needed (F_int acts on >= 4.6e4 rot); old-IC runs are
usable for deep analysis only after rot ~70.

**ck RCE result 09-25 (merged a6eaac93, /viper/ptmp2/jinma/ckrce_0925):** 1-D global-mean ck RCE (nq2, route B, exact
adiabat) has RCB ~7 bar and a COLD deep: T 3/10/100/200 bar = 2235/2662/3862/4267 K vs the old runs' relaxed
3716/4204/5615/6256 K (~1500-1750 K hotter). The old runs' hot deep is MEMORY of the substellar picket-fence IC (thermal
times 1e4-6e5 rot), not an equilibrium: with F_int = sigma Tint^4 imposed, a hot deep adiabat is inconsistent with
Tint (it corresponds to a larger internal flux) and would cool over >= 1e4 rot. Deep entropy is effectively a free
parameter within feasible run times -> USER DECISION pending (cold RCE vs hot vs evolution-model entropy).
Code heating on the profile: median |H|/emission 6e-3, column -8 % (64-cell test column; not checked at nx1 128).

**USER DECISION 09-25: IC deep entropy = option 1, the cold self-consistent 1-D ck RCE (ckrce_nq2_B.txt).** Cross-check running: interior/evolution-model entropy for Teq 2500 K, 0.66 M_J (1.32 R_J at r_in), Tint 538 K (Thorngren) -> /viper/ptmp2/jinma/interior_0925.

**Interior check 09-25 (/viper/ptmp2/jinma/interior_0925):** entropies in k_B/baryon (SCvH zero, MC14-anchored).
RCE start S = 9.86 (3862 K at 100 bar) = cold edge of the interior range; old hot deep S = 12.7 (would need R > 3 RJ).
Literature: UHJ analogues S 9.95-10.62 (3935-4420 K at 100 bar; their M,R quoted from memory, check before citing);
TGF19 Fig. 2 at Teq 2500 & g 942: S 10.42 (4286 K at 100 bar). TGF19 Tint = eps^(1/4) Teq (= our get_Tint, 538 K);
Sarkis+21 Tint 815-879 K at 2500 K. Through our ck atmosphere Tint 538 K gives S ~0.55 below TGF19 (~430 K colder at
100 bar). Our domain radius 1.417 RJ at 1 bar. RECOMMENDATION: keep the RCE start. For a literature deep (S 10.4): do
NOT splice an adiabat (RCB would jump to ~1 bar); rerun ckrce_0925/rce.py with a higher Tint and add a problem/Tint
override (get_Tint hard-wired at deep_hot_jupiter_rt.cpp 1189, 2126, 2343).
