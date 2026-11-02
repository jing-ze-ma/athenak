# NOTE 2026-11-02 viper -> Caltech (USER): heads-up: the dhj / two-stream translation units will be SPLIT on fin51-r-fix (build time 180 s -> < 1 min per edit); do you have unmerged work touching these files?

User asked for fast fix-test builds. Measured: src/pgen/deep_hot_jupiter_rt.cpp (with header-only utils/two_stream_rt.hpp)
was ONE 181 s / 0.9 GB GPU TU (141 s device codegen: 34 TsrtCkChain instantiations + inline-header kernels compiled in
every includer). Branch tusplit-1102 (test, 5e524b3b1 on 0abb64b6), pure code motion, bitwise on GPU + CPU:
- utils/two_stream_rt_pass.cpp, utils/two_stream_rt_chain.hpp + 6 instantiation TUs two_stream_rt_chain_*.cpp,
  utils/two_stream_rt_driver.cpp (picket_fence_two_stream_RT, CkImplStep, CkCad*, ...), utils/two_stream_rt_col3.cpp;
- pgen/deep_hot_jupiter_rt.hpp (private header) + deep_hot_jupiter_rt.cpp (UserProblem) + _hse.cpp + _src.cpp + _diag.cpp;
- 7 pass lambdas -> named functors TsrtX1F..TsrtMIXF (same bodies); new TUs only for PROBLEM dhj / box / red_giant.
Incremental rebuild after an edit: 18-46 s (was ~180 s); clean -j48 286 s (was 351 s).

It is being rebased onto aea58d13b now (tusplit2-test-1102), re-gated, then merged into fin51-r-fix.
QUESTION: do you have unmerged branches editing deep_hot_jupiter_rt.cpp, two_stream_rt.hpp, wbtan_common.hpp,
roche_common.hpp or two_stream_ck_rst*.hpp (e.g. wbtan-rstfix-1101 remainder, mhdspeed rework, r1m/hydro setup)?
If yes, please list them; we either merge them first or I port them onto the split for you. Please NOTE back
before your next pgen edit; I will not merge the split until you answer or ~6 h pass.
