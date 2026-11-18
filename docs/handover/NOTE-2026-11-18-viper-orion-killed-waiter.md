# NOTE 2026-11-18 viper -> Orion (USER): viper killed one of your build waiters by mistake (sorry)

At 04:16 CET (19:16 PST 11-17) a viper agent cleared its own `while pgrep` waiters on the viper login with a pattern
match, and that also killed your waiter PID 2769737 (orion_stag2 build chain). Only the waiter died: the g65jfb build it
waited on was not touched. **The queued g65khe and g65k builds will not start by themselves.** Please re-issue from
/viper/ptmp2/jinma/orion_stag2 (viper did not, to avoid a double build):

    bash -c 'while pgrep -u jinma -f "fastbuild_viper.sh int_g64.* g65jfb" >/dev/null; do sleep 10; done; env J=12 LABEL=g65khe ./fastbuild_viper.sh int_g64he g65kfb.tgz; env J=12 PROBLEM=built_in_pgens LABEL=g65k ./fastbuild_viper.sh int_g64 g65kfb.tgz'

Viper agents now kill only PIDs they started themselves.

-- viper
