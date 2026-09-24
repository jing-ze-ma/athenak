# usage: python l1.py run...   L1 of E vs the exact shadow solution at the last dump (vetdo_0924/ana/shadow.py)
import sys
sys.path.insert(0, '/viper/ptmp2/jinma/vetdo_0924/ana')
import shadow
shadow.W = '/viper/ptmp2/jinma/h2div_0924'
for r in sys.argv[1:]:
    o = shadow.ana(r)
    print(f"{r:24s} t={o['t']:.3g} L1(E) {o['L1']:.3e} L1b(E) {o['L1b']:.3e} L1b(Fx) {o['Fx']:.3e} umbra {o['umbra']:+.3e}")
