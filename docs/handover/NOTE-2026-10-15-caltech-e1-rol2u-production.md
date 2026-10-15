# NOTE Caltech -> Raven (10-14 21:17 PDT): AG Car E1 production is now the ROL2U chain agcE1r (USER); agcE1s stopped

- agcE1r 4553459 (link 1, 12 h backfill slot, then 4553460/61 7-d with time-min 12 h), rol2u_1015/e1/chain, binary 5ab8c87e,
  from agcar3d.00052.rst (t 2.6e7). After 1000 cycles (t 2.6707e7): 0.77 s/cycle (agcE1s ~3.0), dt 727.27 flat, Picard mean 5.0-5.5
  (agcE1s 37-64), 0 DT_RETRY, 0 NON-CONVERGED, 0 FATAL/NaN, 525 hesdirk2 NOT ADMISSIBLE -> BE redos, fofc at-floor 0.
  L_top 1.261e38 (agcE1s 1.260-1.269e38 over 2.6-2.65e7), Mdot_top 2.3e15 (agcE1s 1.8-2.0e15).
- agcE1s 4513330-32 stopped at t ~2.68e7 (rst up to 00053 kept in e1s_1014/prod/rst). The e1 gate DIFF after DT_RETRY (f6473601)
  is still open for you.
