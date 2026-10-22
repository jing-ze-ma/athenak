# NOTE 2026-10-22 caltech -> viper (USER): E2 = accept_fallback + time2_posfix=4 passes the RSG E 129.8 d collapse (re 01599ea3)

Arm E2 (ff1f52bc, md5 cd998ac0): implicit_face_opac=harmonic, time2_vet_col=extrap2, time2_accept_fallback=true, time2_posfix=4;
from prodF rsg.00025 (125 d), 1 node H200 / 4 ranks, debug-QOS segments (rst via rst/rsg.00026 at 128.7 d), tlim 1.1225e7 = 129.92 d.
 dt vs the production (harmonic only) at the same t:
   125.5 d 582 vs 125 | 126.5 d 409 vs 67 | 127.0 d 951 vs 53 | 128.0 d 1000 vs 178 | 128.7 d 400 vs 164 s
   129.8 d (production collapse to 2-5 s, first negative E): E2 min dt 265 s at 129.84 d (one cycle), back to 1000 s; reaches 129.92 d.
 125 -> 129.92 d: 62 DT_RETRY (mostly min(E,T) < 0, each recovered within a cycle), 110 posfix re-solves, 0 fatal; implicit_vimp
 positivity fallbacks 232 in the last segment. No MHD / am_vector in these decks (posfix not wired there).
Not run (cancelled, USER): D2 (retry fix alone) and C2 (extrap2 alone) through the collapse; kirchhoff arms.
Caltech production config (both RSG E from rsg.00021 / 105 d and HG from hegiant.00003 / 15 d), launched in the first allocation that starts:
  ff1f52bc + harmonic + extrap2 + time2_accept_fallback=true + time2_posfix=4.
USER: please consider time2_posfix (N = 4) together with the accept_fallback default flip (01599ea3) -- same gates; and both in the
accretor production if all pass.
