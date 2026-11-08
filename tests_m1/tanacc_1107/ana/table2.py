"""tanfix-1109: per (geometry, level, case) one row per arm, errors vs EXACT (runs2/*/eval.json)."""
import json, os, glob
ROOTS = ['/viper/ptmp2/jinma/openitems_1107/tanacc/runs2', '/raven/ptmp/jinma/openitems_1107/tanacc/runs2']
ARMS = ['off', 'on', 'off_tw', 'on_tw', 'on_twj', 'on_hd', 'C_off', 'C_on', 'C_on_tw']
print('| g lev case | region | arm | E_L1 | E_bias | Fr_L1 | Q_ex | Q_used_err | Drr_L1 | vlat_err | Jfs_L1 |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for g in ('sp', 'cs'):
    for lev in ('lo', 'hi'):
        for case in ('star', 'shd', 'semi', 'mid', 'hot5', 'hot'):
            reg = 'all' if case == 'star' else 'aff'
            for a in ARMS:
                fs = [f'{r}/{g}_{lev}_{case}_{a}/eval.json' for r in ROOTS]
                f = next((x for x in fs if os.path.exists(x)), fs[0])
                if not os.path.exists(f):
                    continue
                d = json.load(open(f)).get(reg)
                if d is None:
                    continue
                print(f"| {g} {lev} {case} | {reg} {d['N']} | {a} | {d['E_L1']:.4f} | {d['E_bias']:+.4f} | "
                      f"{d['Fr_L1']:.4f} | {d['Qex_mean']:.1e} | {d['Q_used_err']:.1e} | {d['Drr_L1']:.1e} | "
                      f"{d['vlat_err']:.1e} | {d['Jfs_L1']:.3f} |")
