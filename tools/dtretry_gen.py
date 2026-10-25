#!/usr/bin/env python3
"""fin4-dtretry: generate the member visitors of the dt-retry checkpoint from clang's AST.

usage: dtretry_gen.py <clang++> <flags file> <src dir>
Writes src/rad_m1/rad_m1_dtretry_visit.hpp, src/hydro/hydro_dtretry_visit.hpp and
src/mhd/mhd_dtretry_visit.hpp, src/diffusion/resistivity_dtretry_visit.hpp (mhdfix-1025): for
every data member of RadiationM1 / Hydro / MHD / Resistivity (and of the nested plain
structs listed in NESTED), one
call of the visitor in declaration order --
  v.view(name, x)    any Kokkos::View (DvceArray*, HostArray*, raw Kokkos::View)
  v.dual(name, x)    any Kokkos::DualView (DualArray*)
  v.value(name, x)   everything copy-assignable that is not a View (scalars, fixed arrays,
                     enums, std::vector/std::string/std::map, plain structs)
  (nothing)          pointers, const members, threads/futures/timers/execution spaces,
                     MPI handles, task ids -- listed in a trailing comment.
Re-run after adding members to these classes; tst/dtretry checks the files are current."""
import json, re, subprocess, sys, os

CLANG, FLAGS, SRC = sys.argv[1], sys.argv[2], sys.argv[3]
SKIP_RE = [r'^std::thread', r'^std::future', r'Timer$', r'^MPI_', r'DevExeSpace', r'TaskIDs$',
           r'^VgdHlPend$', r'^VgdRag', r'^std::vector<MPI_Request>', r'^FofcReport$']
# members that count the exchanges of a peer-to-peer protocol (CUDA IPC halos): the peers'
# flags in device memory have moved on with the failed attempt, so these stay monotonic
# (rolled back, the next exchange would wait for -- or accept -- a stale flag)
# fin51-1014: acc_nccl is the lazily made NCCL communicator's state (with acc_ncomm, a
# pointer, not visited): rolled back to -1 a retry would make a second communicator
SKIP_NAMES = {'hi_seq', 'vgd_ipc_seq', 'acc_nccl'}
NESTED = ['VgdTabHost', 'VgdTabDev', 'M1OpacTab']   # structs holding Views: recurse

def ast(header, name):
    tu = '#include "athena.hpp"\n#include "mesh/mesh.hpp"\n#include "%s"\n' % header
    p = subprocess.run([CLANG] + open(FLAGS).read().split() + ['-x', 'c++', '-fsyntax-only',
                       '-Xclang', '-ast-dump=json', '-Xclang', '-ast-dump-filter=' + name, '-'],
                       input=tu, capture_output=True, text=True)
    s = p.stdout; dec = json.JSONDecoder(); i = 0
    def find(o):
        if isinstance(o, dict):
            if o.get('kind') in ('CXXRecordDecl',) and o.get('name') == name and o.get('completeDefinition'):
                return o
            for v in o.get('inner', []) or []:
                r = find(v)
                if r: return r
        return None
    while i < len(s):
        j = s.find('{', i)
        if j < 0: break
        o, i = dec.raw_decode(s, j)
        r = find(o)
        if r: return [(f['name'], f['type']['qualType'], f['type'].get('desugaredQualType', ''))
                      for f in r.get('inner', []) if f.get('kind') == 'FieldDecl']
    raise SystemExit('no complete definition of %s found (stderr: %s)' % (name, p.stderr[-2000:]))

def isview(t, d):
    return (re.match(r'^(DvceArray|HostArray|ScrArray)\d', t) or t.startswith('Kokkos::View<')
            or 'Kokkos::View<' in d and 'DualView' not in d and not t.startswith('std::'))
def isdual(t, d): return t.startswith('DualArray') or 'DualView' in t
def emit(fields, pre, out, skipped):
    for n, t, d in fields:
        nm = pre + n
        base = re.sub(r'(\[\d+\])+$', '', t).strip()
        dims = re.findall(r'\[(\d+)\]', t[len(base):])
        if n in SKIP_NAMES:
            skipped.append(nm + ' : ' + t + (' (communicator state, kept)' if n == 'acc_nccl' else ' (protocol counter, kept monotonic)')); continue
        if t.startswith('const ') or t.endswith('*') or '(*)' in t or t.endswith('&'):
            skipped.append(nm + ' : ' + t); continue
        if any(re.search(p, base) for p in SKIP_RE):
            skipped.append(nm + ' : ' + t); continue
        if base.startswith('DvceFaceFld'):
            for c in ('x1f', 'x2f', 'x3f'):
                out.append('  v.view("%s.%s", %s.%s);' % (nm, c, nm, c))
            continue
        if base.startswith('DvceEdgeFld'):   # mhdfix-1025: MHD efld
            for c in ('x1e', 'x2e', 'x3e'):
                out.append('  v.view("%s.%s", %s.%s);' % (nm, c, nm, c))
            continue
        if base in NESTED:
            emit(SUB[base], nm + '.', out, skipped); continue
        if base.startswith('std::vector<') and (isview(base[12:-1], '') or base[12:-1].startswith('DvceArray')):
            out.append('  v.value("%s", %s);' % (nm, nm))
            out.append('  for (std::size_t q_ = 0; q_ < %s.size(); ++q_) {v.view("%s[]", %s[q_]);}' % (nm, nm, nm))
            continue
        kind = 'dual' if isdual(base, d) else ('view' if isview(base, d) else 'value')
        if dims and kind != 'value':
            idx = ''.join('[q%d_]' % k for k in range(len(dims)))
            loops = ''.join('for (int q%d_ = 0; q%d_ < %s; ++q%d_) ' % (k, k, dims[k], k) for k in range(len(dims)))
            out.append('  %s{v.%s("%s[]", %s%s);}' % (loops, kind, nm, nm, idx))
        else:
            out.append('  v.%s("%s", %s);' % (kind, nm, nm))

def gen(header, cls, outfile, guard, ns):
    out, skipped = [], []
    emit(ast(header, cls), '', out, skipped)
    txt = ['// GENERATED by tools/dtretry_gen.py from clang\'s AST -- DO NOT EDIT BY HAND',
           '// (fin4-dtretry: the member visitor of the dt-retry checkpoint, driver/dt_retry.cpp)',
           '#ifndef %s' % guard, '#define %s' % guard, '', 'template <class V>',
           ('void %s::%s::DtRetryVisit(V &v) {' % (ns, cls) if ns else
            'void %s::DtRetryVisit(V &v) {' % cls)] + out + ['}', '', '// not visited (pointers, const, threads, timers, MPI, ...):'] + \
          ['//   ' + s for s in skipped] + ['', '#endif  // %s' % guard, '']
    open(os.path.join(SRC, outfile), 'w').write('\n'.join(txt))
    print('%s: %d visits, %d skipped -> %s' % (cls, len(out), len(skipped), outfile))

SUB = {n: ast('rad_m1/rad_m1.hpp', n) for n in NESTED}
gen('rad_m1/rad_m1.hpp', 'RadiationM1', 'rad_m1/rad_m1_dtretry_visit.hpp', 'RAD_M1_RAD_M1_DTRETRY_VISIT_HPP_', 'radm1')
gen('hydro/hydro.hpp', 'Hydro', 'hydro/hydro_dtretry_visit.hpp', 'HYDRO_HYDRO_DTRETRY_VISIT_HPP_', 'hydro')
gen('mhd/mhd.hpp', 'MHD', 'mhd/mhd_dtretry_visit.hpp', 'MHD_MHD_DTRETRY_VISIT_HPP_', 'mhd')
gen('diffusion/resistivity.hpp', 'Resistivity', 'diffusion/resistivity_dtretry_visit.hpp',
    'DIFFUSION_RESISTIVITY_DTRETRY_VISIT_HPP_', '')
gen('hydro/hydro_am.hpp', 'AmModule', 'hydro/hydro_am_dtretry_visit.hpp', 'HYDRO_HYDRO_AM_DTRETRY_VISIT_HPP_', 'hydro')
