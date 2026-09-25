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
