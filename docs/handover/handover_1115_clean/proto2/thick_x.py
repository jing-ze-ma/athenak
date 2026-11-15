"""thick_x.py: extra thick-limit alpha designs (A1^2 with/without reconstruction)."""
import tests as T
T.test_thick([dict(cand='A1sq', recon=None, Dguard=True), dict(cand='A1sq', recon='minmod', Dguard=True),
              dict(cand='A1sq', recon='mc', Dguard=True)])
