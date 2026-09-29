# data-w121-mhd-pkg (data only, not code)

`w121_mhd_etamax_pkg_slim.tgz` (26.7 MB, md5 7f02a2c1b8dd1be4830e2add7920d8fd) = the package for
`docs/handover/TASK-2026-09-29-deltaai-w121-mhd-etamax.md` on rt-integration, WITHOUT the 135 MB reference restart
dhj.00600.rst (not needed: every arm is a fresh start `-i` from remap.dat + the input template).

    git fetch fork data-w121-mhd-pkg
    git show fork/data-w121-mhd-pkg:w121_mhd_etamax_pkg_slim.tgz > w121_mhd_etamax_pkg_slim.tgz
    md5sum w121_mhd_etamax_pkg_slim.tgz   # 7f02a2c1b8dd1be4830e2add7920d8fd
    tar xzf w121_mhd_etamax_pkg_slim.tgz  # -> deltaai_pkg/

Delete this branch on the fork after the transfer.
