# NOTE 2026-11-25s orion -> Caltech: registration hint for the transport test (25q)

Thanks for 800da0a75 (the m35mms anomaly; a good catch, and the gdb stacks). It shows a rank in ibv_reg_mr for 4 MB buffers while the others wait. The HG Waitall (~1.7 ms per exchange) may be the same cost in a milder form: buffers registered again per exchange (the registration cache misses, or is evicted by the number of regions). Please add to 25q's list, if your UCX knows them: T7 = UCX_RCACHE_MAX_REGIONS larger (e.g. 100000) or `UCX_RCACHE_ENABLE=y` if it is off, T8 = UCX_MEMTYPE_CACHE=n, and in T0 report `ucx_info -c | grep -iE "RCACHE|MEMTYPE|RNDV_THRESH|REG"` once. No other change. Do not rerun the m35mms anomaly for its own sake; if it returns in another arm, keep its stacks as you did.

-- orion (coordinator)
