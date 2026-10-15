# NOTE Caltech -> Raven (10-14 20:55 PDT): ROL2U E1 chain agcE1r submitted despite the E1 gate DIFF (USER choice)

- USER chose to run the E1 ROL2U chain anyway after the e1 gate (4552183) showed GATE-BASE / GATE-RETRY DIFF starting at the first
  DT_RETRY (1e-13 at hst row 6 -> 1e-7 hydro, 1.8e-3 m1_f2 after 20 cycles; C = L2 BITWISE; details f6473601).
  Please check the retry restore in 1bd0c2ce (rad_m1_dtretry_visit.hpp suspected): P (5b7e1bc5) vs B (1bd0c2ce default keys) should be
  bitwise through a natural retry.
- agcE1r 4553459 -> 4553460 -> 4553461 (afterany), 1 node x 2 H200, normal QoS, -t 7-00:00:00 --time-min=12:00:00, tlim 1e9,
  binary fin51-r-hdr 5ab8c87e (= 1bd0c2ce + rst header search limit 1 MB, parameter_input.cpp only; md5 0677cd3b), seed agcar3d.00052.rst
  (t 2.6e7). hgR2w 4552184-86 also uses 5ab8c87e. agcE1s production continues; it is stopped only after agcE1r runs cleanly (user plan).
- E1W (narrow wedge) and the He giant LARGE FIN5.1 gate C were stopped by the USER to free GPUs (E1W restartable from
  e1w_1014/prod/rst 00046; large gate C result in 5ae50eb2).
