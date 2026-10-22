# TASK 2026-10-22 caltech -> viper (USER): harmonic FRESH productions hit a dt collapse too (later than with arithmetic)

The fresh starts with implicit_face_opac = harmonic (binary he_star_m1_harm_d3df1401_nofma, = b8f75bfa physics, no extrap2):
- HG global (hgglob_1021/prodF, deck hgglob_fresh_clip_harm.athinput): CFL dt 125 s until t = 1.62e6 s (18.8 d; the arithmetic
  onset was ~5 d). First DT_RETRY cycle 12989 (time2_reject, resid 8.6e-4). By t = 1.7357e6 s (20.1 d) dt is ~1e-4 s and stuck:
  repeated "hesdirk2 stage 1 NOT ADMISSIBLE (converged=1 min(E,T)<0)".
- RSG E (rsg_E_1020/prodF, deck rsg_E_1020_clip_harm.athinput): CFL dt 1000 s until t = 9.26e6 s (107 d; arithmetic ~49 d).
  First DT_RETRY cycle 9260 (time2_reject, resid 1.5e-6). By t = 1.122e7 s (130 d) dt ~2-5 s, same NOT ADMISSIBLE message.
- Summary NON-CONVERGED=[1-9]: 0 in both (rejections are admissibility, not convergence).
USER: both jobs cancelled (5245302, 5240541); chain links 5245303/04, 5240542/43 HELD. No tests on Caltech (USER rule).
Restarts kept: /resnick/groups/carnegie_poc/jingze/hgglob_1021/prodF/rst/hegiant.00004.rst, rsg_E_1020/prodF/rst/rsg.00025.rst
(+ run.log in each prodF).
ASK: diagnose where min(E,T) < 0 appears (photosphere column like the arithmetic case? pole/seam? panel imprint?) and whether
kirchhoff or harmonic + extrap2 (a002299e) avoids it from the last pre-onset restart. The default flip (9c432f17) should wait for this.
