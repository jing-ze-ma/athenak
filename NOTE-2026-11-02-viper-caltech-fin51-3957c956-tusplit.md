# NOTE 2026-11-02 viper -> Caltech (USER): fin51-r-fix = 3957c9562: dhj / two-stream TU split merged (pure code motion, bitwise); rebase dhj/two-stream work onto it

fin51-r-fix aea58d13b -> 3957c9562 = merge of tusplit2-1102 (7992797f9), same tree as 7992797f9. Commits 637bd61f2
(pass + chain/grey launch templates in own TUs), 24531e6e1 (dhj pgen split: deep_hot_jupiter_rt.hpp private header +
deep_hot_jupiter_rt.cpp UserProblem, _hse.cpp, _src.cpp, _diag.cpp), 7992797f9 (heavy inline launchers out of line:
two_stream_rt_driver.cpp, two_stream_rt_col3.cpp). New TUs only for PROBLEM dhj / box_convection / red_giant.

- Rebased by re-running the code-motion scripts on aea58d13b (they reproduce the 0abb64b6 split byte for byte); the
  same 408 distinct changed lines on both bases; all 207 lines of the followup present (rthst DhjFluxHistory hunks now
  in deep_hot_jupiter_rt_diag.cpp, RtFhRstAfterAlloc in two_stream_rt_pass.cpp, flux_hst_rst keys in
  deep_hot_jupiter_rt.cpp).
- Gates vs an aea58d13b binary: GPU fresh r1m 40 cycles and W121 rot-320 restart +20: rst data + mhd.hst IDENTICAL;
  CPU C8 300 cycles and 150+150 (flux_hst_rst on): rst + mhd.hst + user.hst IDENTICAL; box_convection, red_giant build.
- Build: dhj pgen TU 181 s -> slowest TU 43 s; incremental after a dhj/wbtan edit 18-29 s, roche/two_stream_rt 46 s.
  NOTE your incremental build trees: the first build after this merge recompiles the new TUs.

Please: base ck_strang and any dhj-pgen / two_stream_rt edits on 3957c9562. Where code moved: grep the function name
in src/pgen/deep_hot_jupiter_rt*.cpp / .hpp and src/utils/two_stream_rt*.cpp. halofc-1102 is being merged onto
3957c9562 and gated on MI300A now (bvals only, no overlap).
