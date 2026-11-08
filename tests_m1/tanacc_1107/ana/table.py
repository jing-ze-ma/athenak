"""tanacc-1107: markdown tables from runs/*/eval.json (on vs off side by side)."""
import json, glob, os
R = os.environ.get('TA_RUNS', '/viper/ptmp2/jinma/openitems_1107/tanacc/runs')
ev = {os.path.basename(os.path.dirname(f)): json.load(open(f)) for f in glob.glob(R + '/*/eval.json')}
def row(base, reg, keys):
    on, off = ev.get(base.replace('@', 'on')), ev.get(base.replace('@', 'off'))
    if on is None or off is None or reg not in on:
        return None
    cells = []
    for k in keys:
        if k in on[reg]:
            cells.append('%.2e / %.2e' % (on[reg][k], off[reg][k]))
        else:
            cells.append('-')
    return '| %s | %s | %d | %s |' % (base.replace('_@', ''), reg, on[reg]['N'], ' | '.join(cells))
keys = ['E_L1', 'E_bias', 'Fr_L1', 'Flat_L1', 'E_vs_Jfs_L1', 'Jfs_L1', 'Q_used_err', 'Qex_mean', 'Q_fs_err',
        'vlat_err', 'Drr_L1']
print('| run | region | N | ' + ' | '.join(k + ' on/off' for k in keys) + ' |')
print('|' + '---|'*(3 + len(keys)))
for g in ('sp', 'cs'):
    for case in ('star', 'shd', 'semi', 'hot'):
        for lev, sfx in (('lo', ''), ('hi', ''), ('lo', '_hd')):
            base = f'{g}_{lev}_{case}_@{sfx}' + os.environ.get('TA_SFX', '')
            for reg in (('all',) if case == 'star' else ('free', 'aff', 'aff_far')):
                r = row(base, reg, keys)
                if r:
                    print(r)
