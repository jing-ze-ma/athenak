# NOTE 2026-10-22 caltech -> viper (USER): we are testing the harmonic dt collapse HERE (re 951df580)

USER: Caltech runs the collapse test itself; viper need not run its own diagnosis of 951df580 (no duplicate work).
Running: RSG E from rsg.00021 (105 d, first DT_RETRY at 107 d), 2 nodes H200 each, tlim ~112 d:
 A harmonic (d3df1401 binary, control)  B kirchhoff (same binary)
 C ff1f52bc harmonic + time2_vet_col = extrap2   D ff1f52bc harmonic only (isolates prod-retry-1022)
Why D: both production collapses repeatedly take the BE-redo path ("not cured by 2 cuts ... dt_retry_t2_be") with
dt bouncing back up within a cycle, i.e. the path prod-retry-1022 fixes. Also locating the first min(E,T) < 0.
Results -> a NOTE here when done. HG (pre-onset rst hegiant.00003 ~15 d) waits for the RSG verdict.
The default flip (9c432f17) should still wait for this.
