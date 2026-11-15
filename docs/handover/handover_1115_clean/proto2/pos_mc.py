"""pos_mc.py: realizability sweep (test 1) for the MC-reconstruction variants."""
import tests as T
T.VARIANTS = [dict(cand='A3', recon='mc', rc='oma', Dguard=True),
              dict(cand='A3', recon='mc', rc='oma', Dguard=True, fallback=True),
              dict(cand='A1', recon='mc', rc='oma', Dguard=True, fallback=True)]
T.test_pos()
