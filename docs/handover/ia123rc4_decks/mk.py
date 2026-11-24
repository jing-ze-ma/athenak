# mk.py REV (MERGE123): IA12.3-rc3 Caltech decks {tile,notile}/IA123_<star>.athinput + IA123_m35zams_wallsym.athinput (HELD)
#   + IA123_t9l62_restart.args, from merge122/caltech_ia122, in the RELEASE names (<rad_vet>, output variable vet*, problem vet_*),
#   with ONLY the IA12.3 best keys added (COORD rc1): <rad_vet> vet_gd_wall_sym, vet_gd_theta_seam_sym, clean_clos1_budget,
#   clean_bud_maxnorm = true, clean_vel_settled = false (named: breaks wb_sp / wb_cs), <hydro> floor_from_ekin = true.
#   (rad_signal_speed_diff and max_stage2 are already in the IA12.2 decks.)
import sys, os, re, hashlib, difflib
REV = sys.argv[1]; B = '/orion/ptmp/jinma/clean_1115/stag2/merge122/caltech_ia122'
RV = ['vet_gd_wall_sym = true   # IA12.3 WALLSYM: mirror-paired phi-wall maps (identical on whole-period wedges)',
      'vet_gd_theta_seam_sym = true   # IA12.3: theta seam face from both sides (equatorial mirror); FATAL without it on theta-periodic meshes',
      'clean_clos1_budget = true   # IA12.3-rc3 HALFDAY-M35 L1 + RWFIX floor (radwave 2.4-5.5e-5): stage-1 pass 0 in the budget norm',
      'clean_bud_maxnorm = true   # IA12.3-rc3 HDM35-C 83b71abed form (HALF12): budget test in the per-cell max norm',
      'clean_vel_settled = false   # named OFF: HALF12 drops it (the correct form costs more than it saves)']
HY = ['floor_from_ekin = true   # IA12.3 FLOORT/GEOMAUDIT/LZCUT: c2p keeps E; vceil KE -> heat; KE-dominated tfloor deficit paid from own r/theta KE (phi kept, L_z exact)']
GONE = ['implicit_pos_floor_solve', 'implicit_fast_kernels', 'implicit_allow_multid', 'cs_cross_mono', 'vet_gd_reach', 'clean_nl_print']
def rename(L):
    O = []; blk = ''
    for l in L:
        m = re.match(r'\s*<(\w+)>', l)
        if m:
            blk = m.group(1)
            if blk == 'rad_m1': l = l.replace('<rad_m1>', '<rad_vet>'); blk = 'rad_vet'
        elif blk == 'problem':
            l = re.sub(r'^(\s*)m1_', r'\1vet_', l); l = re.sub(r'^(\s*pgen_name\s*=\s*)rad_m1_', r'\1rad_vet_', l)
            l = re.sub(r'^(\s*)he_ic_linhao(_|\s*=)', r'\1he_ic_profile\2', l)   # (COMPAT: the new names of the bsg IC keys)
        elif blk == 'rad_vet' and re.match(r'\s*(%s)\s*=' % '|'.join(GONE), l):
            continue   # (COMPAT KL5: "no effect" on the trunk; stripped so the start-up WARNING stays clean)
        elif blk.startswith('output'):
            l = re.sub(r'^(\s*variable\s*=\s*)m1(\b|_)', r'\1vet\2', l)
        O.append(l)
    return O
def add(L, tag, rv=RV, hy=HY):
    O = []; blk = ''; nrv = nhy = 0
    have = {k.split()[0] for k in rv + hy if any(re.match(r'\s*%s\s*=' % k.split()[0], l) for l in L)}
    for l in L:
        m = re.match(r'\s*<(\w+)>', l)
        if m: blk = m.group(1)
        if any(re.match(r'\s*%s\s*=' % k, l) for k in have): continue   # replaced below
        O.append(l)
        if blk == 'rad_vet' and re.match(r'\s*vet_gd_wall_interp\s*=\s*true', l): O += rv; nrv += 1
        if m and blk == 'hydro' and hy: O += hy; nhy += 1
    assert nrv == 1 and nhy == (1 if hy else 0), (tag, nrv, nhy); return O
for v in ['tile', 'notile']:
    os.makedirs(v, exist_ok=True); md = []
    for s in ['bsg', 'm35mms', 't9l62']:
        L = open(f'{B}/{v}/IA122_{s}.athinput').read().splitlines()
        O = [f'# IA12.3-rc3 2026-11-24 ia123-rc3 {REV} (release trunk vet-clean-release + IA12.3). = caltech_ia122/{v}/IA122_{s}.athinput in the release names (<rad_vet>, vet*) with ONLY the IA12.3 best keys (mk.py header)'] + add(rename(L), s)
        f = f'{v}/IA123_{s}.athinput'; open(f, 'w').write('\n'.join(O) + '\n')
        open(f'{v}/diff_vs_IA122_{s}.txt', 'w').writelines(difflib.unified_diff([x + '\n' for x in L], [x + '\n' for x in O], f'IA122_{s}', f'IA123_{s}', n=0))
        md.append(f'{hashlib.md5(open(f, "rb").read()).hexdigest()}  IA123_{s}.athinput\n')
    open(f'{v}/MD5SUMS', 'w').writelines(md)
L = open(f'{B}/IA122_m35zams_wallsym.athinput').read().splitlines()
O = [f'# IA12.3-rc3 M35ZAMS Caltech deck (HOLD for the USER), {REV}. = caltech_ia122/IA122_m35zams_wallsym.athinput in the release names + the IA12.3 best keys'] + add(rename(L), 'm35z')
open('IA123_m35zams_wallsym.athinput', 'w').write('\n'.join(O) + '\n')
A = open(f'{B}/IA122_t9l62_restart.args').read()
A = A.replace('IA12.2 restart', 'IA12.3-rc3 restart').replace('clean-stag2-interim-f 465aff485', f'ia123-rc3 {REV}').replace('caltech-ia122-465aff485', f'caltech-ia123rc3-{REV}').replace('IA121_t9l62_restart.args', 'IA123_t9l62_restart.args')
A = re.sub(r'(?m)^rad_m1/', 'rad_vet/', A)
A += ('# IA12.3-rc3 keys (release names; the IA11 rst stores <rad_m1>, read as <rad_vet> by the aliases; none of these is stored there):\n'
      + ''.join(f'rad_vet/{k.split()[0]}={k.split()[2]}\n' for k in RV) + ''.join(f'hydro/{k.split()[0]}={k.split()[2]}\n' for k in HY))
open('IA123_t9l62_restart.args', 'w').write(A)
