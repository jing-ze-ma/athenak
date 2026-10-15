# NOTE Raven -> Caltech (10-15 night): ROL2U safety audit -> SAFE WITH CONDITIONS; please add these monitors to agcE1r / hgR2w reports

Raven audit (rol2u_safety_1015): accepted-solve error at the ROL2U stop median 1e-7, p99 3e-5 (E1 stiff states), <= 4e-9 (He giant);
no masked failures; floor energy from the opacity lag ~1e-5 L_top, flat; local zone profiles at the chaos level; retries in BASE are
false alarms (stalls at resid ~5e-8), ROL2U takes 0. Long-run statistics (to 5e8 s) still pending on Raven/Viper.
CONDITIONS (already true for agcE1r/hgR2w -- please keep): use the FULL ROL2U bundle (implicit_tol_rel 0.1 + _max 1e-4 + opac_lag pred/2
passes + gas_newton_u4). NEVER implicit_tol_rel alone (audit: unsafe on post-ramp E1, accepts solves BASE rejects). cs2_vimp_npass 3.
MONITORS per report (vs agcE1s / hgF4w at matching times): "mean tol_eff" and "lagged solves" lines; e_floor clip count + floor energy
created (event log m1_pos_un) vs L_top; NON-CONVERGED and DT_RETRY per 1000 cycles; NOT ADMISSIBLE rate; Picard mean; s/cycle;
E1: d105 / vr105 / rsurf (agcE1r is in an infall-like phase). Flag at once: Picard mean > 15, NC > 1 per 1000 cycles, floor energy > 1e-4 L_top.
