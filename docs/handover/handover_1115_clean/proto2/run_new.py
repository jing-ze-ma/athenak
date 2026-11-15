"""run_new.py: coordinator additions - linear sensor variant S1 and van Leer, side by side.
Usage: python3 run_new.py pos | posrf | thin <name> <chi> | quick | quickrf  [--slin]"""
import sys
import tests as T
V = {'S1': dict(cand='S1', recon=None, npic=1, Dguard=True),
     'S1fb': dict(cand='S1', recon=None, npic=1, Dguard=True, fallback=True),
     'S1raw': dict(cand='S1', recon=None, npic=1),
     'A3vl': dict(cand='A3', recon='vl', Dguard=True),
     'A3vlfb': dict(cand='A3', recon='vl', Dguard=True, fallback=True),
     'A3mc': dict(cand='A3', recon='mc', Dguard=True),
     'A3mcfb': dict(cand='A3', recon='mc', Dguard=True, fallback=True),
     'A3mm': dict(cand='A3', recon='minmod', Dguard=True),
     'A3mmfb': dict(cand='A3', recon='minmod', Dguard=True, fallback=True),
     'A31st': dict(cand='A3', recon=None, Dguard=True),
     'A31stfb': dict(cand='A3', recon=None, Dguard=True, fallback=True),
     'A1mc': dict(cand='A1', recon='mc', Dguard=True),
     'RF': dict(cand='rayflux', recon=None, npic=1, Dguard=True),
     'RFfb': dict(cand='rayflux', recon=None, npic=1, Dguard=True, fallback=True),
     'RFraw': dict(cand='rayflux', recon=None, npic=1)}
import scheme as S
if '--slin' in sys.argv:          # piecewise-linear ray source (needed by rayflux)
    S.SLIN = True
    sys.argv.remove('--slin')
if sys.argv[1] == 'pos':
    T.CFLS = (0.5, 4.0, 100.0, 1000.0)
    T.VARIANTS = [V[k] for k in ('S1raw', 'S1', 'S1fb', 'A3vl', 'A3vlfb', 'A3mc', 'A3mcfb',
                                 'A3mm', 'A3mmfb', 'A31st', 'A31stfb')]
    T.test_pos()
elif sys.argv[1] == 'posrf':
    T.CFLS = (0.5, 4.0, 100.0, 1000.0)
    T.VARIANTS = [V[k] for k in ('RFraw', 'RF', 'RFfb')]
    T.test_pos()
elif sys.argv[1] == 'quickrf':
    sel = [V['RF']]
    T.test_thick(sel)
    T.test_oddeven(sel)
    T.test_cons(sel)
    T.test_torder_thick(sel)
    T.test_torder(sel)
    T.test_rwave(sel)
elif sys.argv[1] == 'thin':
    T.test_thin([V[sys.argv[2]]], chis=(float(sys.argv[3]),))
elif sys.argv[1] == 'quick':
    sel = [V['S1'], V['A3vl'], V['A3mc']]
    T.test_compact
    T.test_thick(sel)
    T.test_oddeven(sel)
    T.test_cons([V['S1'], V['A3vl']])
    T.test_torder_thick([V['S1'], V['A3vl']])
    T.test_torder([V['S1'], V['A3vl'], V['A3mc']])
    T.test_rwave([V['S1'], V['A3vl']])
