"""FIN3 page (worker PG, 2026-10-13): rewrite idort_results.json for the FIN3 final scheme (fin3-1012 281dd686).
Base = idort_results.json.bak_fin3 (the published FIN / twinfix state).  For every test side the old 'final' becomes a
short 'before' line and 'final' holds the FIN3 number + verdict, or PENDING when the final-config rerun has not landed.
Numbers: IN's evaluation (fin3_1012/ana/out_281dd686), this page's figure captions (figsrc/plot.py on the FIN3 runs),
RT's NOTE-retard (bubbles reference).  Edit the tables below as runs land and rerun:  python3 update_fin3.py"""
import copy
import json
import os

P = '/raven/ptmp/jinma/idort_page'
B = json.load(open(P + '/idort_results.json.bak_fin3'))
R = copy.deepcopy(B)
BIN = 'fin3-1012 281dd686 (FIN3 final)'
PEND = 'pending (final rerun running)'


def short(x, n=230):
    if isinstance(x, dict):
        v, r = x.get('verdict'), x.get('result')
    else:
        v, r = None, x
    r = '' if r is None else str(r)
    if len(r) > n:
        r = r[:n].rsplit(' ', 1)[0] + ' ...'
    return {'verdict': (v or '').split(' (')[0] or None, 'result': r}


# ---------------------------------------------------------------- per test: order -> side -> (verdict, result)
# verdict 'PENDING' = the final-config run is still running (Caltech / DeltaAI / Viper chains or Raven queue)
T = {
 1: {'cartesian': ('PASS', 'row 1.006 / peak 0.417 / floor 4e-19; Picard 33.3, 6 NC (Raven 281dd686 build; Viper build 29.9, 5 NC; Picard-heavy at c dt/dx 256 with ring 4, as before)'),
     'spherical': ('PARKED', 'T8 parked, not run')},
 2: {'cartesian': ('PASS', '1.006 / 0.417 at r12 / c0.5 / c32; Picard 3.13 / 2.94 / 3.85; cbh: row 1.006, peak 0.417, Picard 4.0')},
 3: {'cartesian': ('PASS', 'row 0.998, peak 0.737, FWHM 5.6 cells; Picard 3.87')},
 4: {'cartesian': ('PASS', 'completes: row 0.89 / 0.82 (y = -1.9 / -1.5 ... 1.9), peak 0.071; Picard 34.1, 6 NC (FIN2 control: same)')},
 5: {'cartesian': ('PASS', 'Gauss 8x32 (256^2): E in 0.06 %, E out 1.5 %, F out 1.1 %, rF/rF_ex 1.009, f_xx(0.5) 0.3331 / 0.3331; 512^2: E in 0.05 %, E out 1.6 %, F out 1.2 %, rF/rF_ex 1.007, f_xx(0.5) 0.3331 / 0.3331 (MI300A plain build) (FIN2 control 512^2: 1.4 % / 1.3 % / 0.993)'),
     'spherical': ('PASS', 'T6b r<0.9 +0.58 / -0.11 / -0.06 %, r>3 +6.4 / -2.2 / -1.2 % at N 64/128/256 (converging), F r^2 0.987 at 256; T6a E L1 3.9e-6, F r^2 1.3e-3')},
 6: {'cartesian': ('PASS', 'L1 1.11 % (eps 1e-1) ... 0.01 % (1e-8), E_top/ex 0.9994-1.0000; Picard 1.4-2.2'),
     'spherical': ('PASS', 'nlte2 r0/w 1000 vs exact -0.25 / +0.003 / +0.08 / +0.02 % per tau band; E_top sp/cart 0.994')},
 7: {'cartesian': ('PASS', 'ring 12x2 / 36x2: L1 0.011 / 0.011, wedge 158.2 (exact 157.6), umbra 0.21 / 0.19; Picard 3.0 / 3.0'),
     'spherical': ('PASS', 'twin_state trans in the scheme: lit L1 T7L1 at r 4.05 / 4.99 / 5.77 (120 dirs; 432 dirs T7L4), umbra T7UM, penumbra L1 0.18-0.19 (120); Picard 2.25')},
 8: {'cartesian': ('PASS', 'L1 0.030, umbra 0.21; Picard 3.0')},
 9: {'cartesian': ('PASS', 'L1 0.02 / 0.06 / 0.28 %, peak 1.0001'),
     'spherical': ('PASS', 'L1 7.2e-3 / 1.8e-3 / 4.6e-4 (X 2.01 / 2.00), T 2.00 / 2.01 / 1.94')},
 10: {'cartesian': ('PASS', 'ddf / ddiff at t = 16: peak 0.9700 (within 5 %), centroid x_c - v t = -1e-4, L1 0.0237, D_eff/D 0.984; identical to the FIN2 control to 4 digits')},
 11: {'cartesian': ('PASS', 'ddn at t = 16: peak 0.9983, centroid x_c - v t = 0.0000, L1 0.0065, D_eff/D 1.001; identical to the FIN2 control')},
 12: {'cartesian': ('EXCLUDED', 'EXCLUDED from the test suite (user decision, 2026-10-13; not counted as PASS / OPEN / PENDING; numbers for information only). Outflow (IDORT BCs), ring 200, vs the open retarded reference at 0.05 / 0.25 / 1 t_c: L1 BUB200, '
                     'E_tot BUBET; staircase reference (OD): ring 100 0.214 / 0.249 / 0.135, ring 200 0.213 / 0.248 / 0.134, Gauss 8x32 k3 0.190 / 0.263 / 0.174; the periodic deck and bubu variants are not evaluated further (excluded). RT (retard-1012): the time-dependent '
                     '(retarded) formal-solution reference reaches 0.213 / 0.214 / 0.160, i.e. the remaining error at 0.25 t_c is light '
                     'that has not arrived yet (steady VET tensor). IDORT\'s own snapshot scores 0.27 at 1 t_c.')},
 13: {'cartesian': ('PASS', 'E 96.79505, T 3.13663 (target 96.795, 3.13663); Picard 1.06'),
      'spherical': ('PASS', 'relax1 / relax2 early: 1e-12 by dt/4 (bath floor)')},
 14: {'cartesian': ('PASS', 'E 145.78780, T 3.47480 (target 145.788, 3.4747); Picard 1.09'),
      'spherical': ('PASS', 'relax2 early: 1e-12 by dt/4 (bath floor)')},
 15: {'cartesian': ('PASS', 'E 1.132601, v 2.955375, P_tot 3.000000; Picard 5.7')},
 16: {'cartesian': ('PASS', 'self-convergence ladders rw_t10a / rw_t1000a: X 1.94-2.07, T 1.6-2.1, C 1.9-2.1; convergence vs the exact eigenmode (tau 0.01 / 1 / 10): RW16'),
      'spherical': ('PASS', 'rsw X dens 1.83, v 1.99, E 2.04, F1 1.99; T 2.00-2.09 every variable')},
 17: {'cartesian': ('PASS', 'self-convergence ladders rw_p10t100 / rw_p10t1: X 1.97-2.25, C 1.91-2.10, T 1.97-2.11 (p10t100) and 1.33-1.98 (p10t1, OPEN ladder row); convergence vs the model root (tau 0.01 / 1 / 100): RW17'),
      'spherical': ('OPEN', 'rsw T 2.0; rswv X TEST-fixed by the user-approved rswv deck change (implicit_f0_init zero + outflow gas top): 1.73-1.95 over N 64-1024, interior 2.0; rswv T OPEN: dens / eint 1.6-2.06, v / E / F1 1.2-1.3 at the 16-32 pair (the per-step gas-radiation O(dt) term of rw_p10t1)')},
 18: {'cartesian': ('PASS', 'L1 rho / T_gas / T_rad vs the IDORT transport reference: M 1.2 (512) 0.0004 / 0.0004 / 0.0003; M 2 (512) 0.0007 / 0.0011 / 0.0014, (1024) 0.0004 / 0.0006 / 0.0007; M 3 (512) 0.0012 / 0.0009 / 0.0008, (1024) 0.0010 / 0.0017 / 0.0017 (all <= 0.2 %, criterion 2 %); identical to the FIN2 control')},
 19: {'cartesian': ('PASS', 'reproduced: IDORT beam (540 a T0^4 per beam): halo -0.208 / -0.280 / -0.392, dense cloud (rho > 3) from -0.065 / -0.034 / 0.011, peak rho 36 / 40 / 33 at t = 0.06 / 0.10 / 0.16; IDORT snapshots -0.205 / -0.275 / -0.390, from -0.065 / -0.035 / 0.015; Picard 3.9. Run: IN\'s FIN3 cloud deck + args with the IDORT beam, beam_e 5.4e-5 (idort_page/reruns/fin3/cloud_FIN3_b540).')},
 20: {'cartesian': ('PASS', 'L1(E) MARSHC1 / 0.0141 (c dt/dx 1 / 10)'),
      'spherical': ('OPEN', 'X dens 1.09, E 0.56 (front, 1st order); T orders 0.2-0.4 (flat, floor)')},
 21: {'cartesian': ('PASS', "probe T at t = 10, c dt/dx 220: 0.938 / 0.856 / 0.772 / 0.686 / 0.550; c dt/dx 2200: 0.939 / 0.858 / 0.774 / 0.688 / 0.554 (dt self-consistency within 0.7 %; FIN2 control 0.4 %).. Caveat (IN): the probe-5 arrival time t(T >= 0.5) is 4.20 (c dt/dx 220) vs 3.31 (2200), against FIN2's 3.81 / 3.71, so FIN3 is less dt-consistent than FIN2 here; T crosses 0.5 slowly there (dT/dt ~ 0.017) and pointwise T differs by at most 2.7 %. PASS is on the end state, and the criterion is our own (no quantitative reference)")},
}
SX = {  # spherical_extra, by test_id prefix
 'P1 patm': ('PASS', 'E/J-1 at R_ph -0.020 / -0.058 / -0.030 %; max|Fr^2-1| 3.8e-4; band 10-20 order 2.07 / 2.32'),
 'P1T patmT': ('PASS', '+0.021 / -0.059 / -0.029 %; X dens 1.98, E 1.86, F1 2.00; Picard 1.3-1.9'),
 'C patmT eps': ('PASS', '-0.017 / -0.029 % (eps 1e-2), -0.029 % (1e-4)'),
 'P2 pinc': ('PASS', 'new definition (wide channel pinc_ha20_sr30, hole 0.2 rad, screen r 2-3) vs the staircase Monte Carlo L_ex: '
                     'G2 480 dirs 1.010 / 1.014 (both direction sets, rule within 5 %); G1 120 dirs with vet_rotate_avg 4 (FIN3) '
                     '0.999 / 1.008; G1 without rotation 0.996 / 0.898 (ray-effect spread). Narrow channel (stress case): 1.03 / 1.24.'),
 'L pulse': ('PASS', 'L1 8.8e-2 (r0/w 1, true focusing) / 3.7e-4 (r0/w 1000); E_peak sp/cart 1.0009'),
 'L rsw': ('PASS', 'L1 1.4e-4 / 2.3e-4 (rsw / rswt, N 128)'),
}
EX = {  # extra_tests, by test_id prefix
 'rw_t10a X': ('PASS', 'X E 1.96, F1 2.06, rho 2.04, v 2.07; C 1.91-2.06 (fin3 ana/out_281dd686 lad.txt)'),
 'rw_t10a T': ('PASS', 'E 2.67 2.00 1.94; rho 1.98 1.99 1.96; v 1.91 1.91 1.76; F1 1.89 1.87 1.62 (l 1-16)'),
 'rw_t1000a': ('PASS', 'X 1.94, T 1.99-2.09, C 1.96-1.99'),
 'rw_p10t100': ('PASS', 'X 2.03, T 1.97-2.11 every variable, C 2.02-2.04'),
 'rw_p10t1': ('OPEN', 'X 1.97-2.25, C 1.91-2.10; T E 1.82 1.92 1.98, rho 1.93 1.86 1.74, v 1.90 1.84 1.73, F1 2.58 1.98 1.33 '
                      '(improved from rho 1.58, F1 1.26; small errors, rho 1.6e-5; within 0.05 of the FIN2 control)'),
 'grey atm': ('PASS', 'E 2.10 / 2.12 / 2.08 (32-512, steady hesdirk2); tau bands as FIN2 (unchanged)'),
 'pulse k': ('PASS', 'X E 2.15 / 2.27 / 1.91, F1 2.15 / 2.60 / 1.90 (k 0.128 / 12.8 / 1280); C E 2.15 / 2.27 / 1.91, k 0.128 C F1 2.14 (k 0.128 / 12.8 / 1280; k 0.128 C finest pair 256-512)'),
 'G5 2-D': ('PASS', 'd(sigma^2)/dt / 2D = 0.99957 / 0.99987 / 1.00000 / 1.00000 (kappa 12.8 / 128 / 1280 / 12800, target 1 +- 0.02); same as FIN'),
 'T9 rsw': ('PASS', 'X dens 1.83, v 1.99, E 2.04, F1 1.99; T 2.00-2.09'),
 'T9v rswv': ('OPEN', 'Deck fix (user-approved, rswv decks only; code defaults unchanged): comoving radiation IC (implicit_f0_init = zero) + outflow gas top (ox1_bc = outflow). X TEST-fixed: full domain dens 1.73 1.95 1.80, v 1.87 1.90 1.84, E 1.89 1.88 1.81, F1 1.87 1.90 1.85 (N 64-1024; FIN3 decks 0.13-0.70, wall layer); r < 2.5 all 2.0. T: dens 2.03 2.06 1.94 1.69, eint 2.02 2.03 1.90 1.59 (FIN3 1.27 1.17 1.13 1.35); v / E / F1 still 1.94 1.82 1.49 1.27 / 1.94 1.82 1.39 1.21 (levels 1-32): the per-step gas-radiation O(dt) term of rw_p10t1, OPEN as a method item.'),
 'T9t rswt': ('PASS', 'X dens 1.95, v 1.76, E 1.77, F1 1.85; T dens 2.56 / 2.27, v 1.82 / 1.94, E 1.36 / 2.20, F1 1.75 / 2.15'),
 'T9tv rswtv': ('OPEN', 'X dens 1.35, v 0.75, E 0.68, F1 0.66; T 1.4-2.3'),
 'rswt_q': ('PASS', '6.2e-5 / 1.5e-5 / 3.8e-6 (2nd order), rho 7.3e-5 -> 6.5e-6'),
 'rlat lateral': ('PASS', 'X 1.89-2.14 (rlat), 1.93-2.10 (rlat_om5); T rlat_T 1.67 / 1.73, rlat_om5_T 1.72 / 1.74 (rotated, start-up)'),
 'rlatt / rlatt_om5': ('OPEN', 'rlatt_om5_T OPEN (user decision): E and F1 count as converged on the finer pairs (criterion, user-approved, levels 4-16): E 1.29 2.04 2.46, F1 1.37 2.05 2.30, F2 1.33 2.05 2.02 (the 1.3 was the coarsest pair); but v_phi 2.20 1.61 1.49 (1.61 / 1.49 on the finer pairs) and dens 2.02 2.03 1.73 (1.73 at the finest pair) stay below 2. OT is checking whether v_phi is the same per-step conservative2 velocity term (as rw_p10t1). rlatt X stays OPEN: E 2.07 1.96, F3 1.93 1.83, dens 2.44 2.03 (32-256), but F1 / v_r sit on a lateral-resolution-independent floor (order 0.77 then 0.0 at 128-256; about 1.5 % of the small mean |F1|, about 2e-7 absolute); not the rotation (off: worse), interp / cut_thin / symref, and only 30 % lower with 18x24 directions.'),
 'L_z conservation': ('PASS', '+4.8e-6 / -3.9e-5 of L_rad; dE_tot -1.9e-4 / -7.9e-5'),
 'T1 relax1_T': ('PASS', 'errors flat ~1e-7 (bath floor)'),
 'T4 dif T': ('PASS', 'X 2.01 / 2.00; T 2.00 2.01 1.94'),
 'T5 Marshak T': ('OPEN', 'X dens 1.09, E 0.56; T orders 0.2-0.4 (flat)'),
 'T6a thin sphere': ('PASS', 'E L1 3.9e-6 at every N (floor); F r^2 1.9e-2 -> 1.3e-3'),
 'slab:': ('PASS', 'front cell E/E_ex 0.988 / 0.998 / 1.003 / 1.000 / 1.000 (N 32-512); gap L1 1.3e-2 -> 1.1e-4; Picard 2.05-2.12'),
 'cold opaque spherical shell': ('PASS', 'cavity first / mid / front cell E/E_free 0.987 / 0.997 / 0.999 (N 32), 0.996 / 0.999 / 1.000 (64), '
                                         '0.999 / 1.000 / 1.000 (128, 256): fixed in FIN3 (was 211 / 109 / 55 / 28)'),
 'single opaque disc': ('OPEN', 'TEST criterion change (user approved): Reference = the exact steady solution of the STAIRCASE disc the grid holds (cells with centre r <= R, black; ray cast / convex-hull projection), which is brighter than the smooth disc (hull perimeter / 2 pi R 1.132 / 1.054 / 1.020 at n256 / 512 / 1024). E outside / exact(STAIR), rho kappa 1e5: 0.784 / 0.873 / 0.927 at n256 / 512 / 1024 (first order in dx/diameter; steady from t 0.005 to 0.2). n512, rho kappa 1e3 / 1e5 / 1e7: 0.820 / 0.872 / 0.873. Two pieces: (1) vet_sc formal-solution geometry (a steep ray steps one x1 column and skips rows: J_FS/J_stair 0.86 in the x2-near sector at n512, 0.98 in the x1-near), (2) the M1 moment step at the staircase corners (E/J_FS 0.94-0.96, 0.92 at rho kappa 1e3). Candidate, NOT adopted: vet_path = cells (open-od-1013 06f87b40) subdivides steep steps at material edges: FS isotropic (0.98 in every sector), E/exact 0.895 / 0.942 / 0.943 at n512 (1e3 / 1e5 / 1e7). The M1 corner deficit (piece 2) stays OPEN.'),
}
EXNOTE = {'single opaque disc': "Reference changed (user approved, 2026-10-13): smooth disc -> staircase exact (the test's own discrete geometry). On the old smooth reference FIN3 read 0.92 (1e5); the staircase is ~5 % brighter at n512, so the number drops to 0.87 although the solution is unchanged."}  # OD patch (open-1013)
EXBIN = {'T9v rswv': 'fin3-1012 281dd686 (FIN3 final), rswv decks: implicit_f0_init = zero, ox1_bc = outflow'}  # OT patch (open-1013): user-approved rswv deck change
NOTES = {
 1: 'Ring 4 at c dt/dx 256 is Picard-heavy (hc thin correction on vacuum beams); ring 12 and c dt/dx <= 32 converge in 3-4. Spherical T8 parked (vet_gd_beam not flux-conservative); the production channel beam is P2.',
 4: 'Completes with the jse face reconstruction; Picard 34 / 6 NC is the known hc vacuum-beam convergence weakness.',
 5: 'T6b remaining error is radial (M1 core luminosity at the opaque -> vacuum surface), converging about first order.',
 6: 'Spherical uses vet_col_source = gas for held-T scattering.',
 7: 'Spherical: the vet_gd twin uses the mean-transmission state (twin_state trans) and vet_rotate_avg 4; lit cells are within 1 % at 120 directions, the remaining error is penumbra resolution.',
 12: ('IDORT outflow radiation BCs (implicit_bc_x2/x3 marshak; vet_sc rotates with the shared clock). Retardation-limited (RT, '
      'retard-1012): the steady VET tensor carries light that has not arrived; the time-dependent reference reaches L1 0.213 / '
      '0.214 / 0.160. EXCLUDED from the suite by user decision (previously PASS on the outflow criterion, L1 <= IDORT-level 0.27 at every time; information only). Opaque discs under-emit (single '
      'disc 0.87 of the staircase exact at n512, OPEN, listed separately). ' + 'Reference: the 20 sources are 1.5-15 cells in radius, so the grid holds staircase discs whose hull width is 1.03-1.22x the disc diameter (+4-5 % luminosity-weighted). With the STAIRCASE open retarded reference (same construction as the single disc) FIN3 reads L1 0.214 / 0.249 / 0.135 (ring 100), 0.213 / 0.248 / 0.134 (ring 200), 0.190 / 0.263 / 0.174 (Gauss 8x32 k3) at 0.05 / 0.25 / 1 t_c, E_tot 0.88 / 0.86 / 0.87 (ring 200); (would be PASS on L1 <= 0.27; test excluded). The smooth-reference numbers (0.208 / 0.248 / 0.114) are kept for comparison. Candidate vet_path = cells (not adopted), ring 100: 0.206 / 0.251 / 0.099, E_tot 0.90 / 0.89 / 0.90 (staircase reference).'),
 15: 'conservative2 conserves P_tot exactly.',
 19: 'IDORT beam (540 a T0^4 per beam); plotted in IDORT\'s convention (jet, log rho 1-40).',
 20: 'Spherical: the front converges at first order (expected); the T ladder is flat at the front floor, as in the Cartesian twin.',
 21: 'Reproduced qualitatively (no quantitative Cartesian reference).',
 22: 'Not feasible (NOTE-suite §3.6).', 23: 'Not feasible (NOTE-suite §3.6).', 24: 'Not feasible (NOTE-suite §3.6).',
}
NEW_EX = [
 {'test_id': '2-source retardation test (Cartesian, outflow; RT)', 'geometry': 'Cartesian', 'X': None, 'T': None,
  'what': 'a bright and a dim source; before the bright light arrives the dim bubble must stay round',
  'final': {'verdict': 'EXCLUDED', 'binary': BIN,
            'result': 'EXCLUDED (user decision: the light-travel / retardation item is deferred; information only). Retardation-limited: dim bubble ellipticity 0.67 / 0.78 / 0.83 / 0.86 at t 0.01-0.04 (exact 0), dim-field '
                      'L1 1.06-1.88 before arrival (dim-only run 0.09-0.21); total L1 0.16 -> 0.05, E_tot 0.95-0.97. RT: a '
                      'time-dependent formal solution reduces the dim-field error to 0.29-0.46.'},
  'notes': 'new entry (RT, retard-1012); EXCLUDED with the multi-source bubbles (user decision); no cheap fix in FIN3 (vet_causal not included)'},
]


def setf(side, v, r):
    old = side.get('final')
    side['before'] = short(old if old else {'verdict': side.get('verdict'), 'result': side.get('result')})
    side['final'] = {'verdict': v, 'result': r, 'binary': BIN}
    side.pop('latest_fix', None)


def fill(R, sub):
    """sub: placeholder -> value (numbers from the figure captions)"""
    for t in R['tests']:
        for s in ('cartesian', 'spherical'):
            d = t.get(s)
            if isinstance(d, dict) and isinstance(d.get('final'), dict):
                for k, v in sub.items():
                    d['final']['result'] = d['final']['result'].replace(k, v)


def main(sub):
    for t in R['tests']:
        o = t['order']
        for s in ('cartesian', 'spherical'):
            d = t.get(s)
            if not isinstance(d, dict):
                continue
            if o in T and s in T[o]:
                setf(d, *T[o][s])
            elif o >= 22:
                setf(d, 'PARKED', 'not run (parked)')
        vs = [((t.get(s) or {}).get('final') or {}).get('verdict') for s in ('cartesian', 'spherical') if isinstance(t.get(s), dict)]
        t['final'] = {'verdict': vs[0] if vs else None, 'binary': BIN}
        t['notes'] = NOTES.get(o, '')
    for t in R['spherical_extra']:
        if t['test_id'].startswith('P2 pinc'):
            t['notes'] = ('P2 redefined (user approved) to the resolved wide channel; the narrow channel (0.03 sr vs ~0.1 sr per '
                          'direction) is kept as an angular stress case. vet_rotate_avg 4 removes the 120-direction set dependence.')
        for k, (v, r) in SX.items():
            if t['test_id'].startswith(k):
                t['before'] = short(t.get('final'))
                t['final'] = {'verdict': v, 'result': r, 'binary': BIN}
                t.pop('latest_fix', None)
    for t in R['extra_tests']:
        for k, (v, r) in EX.items():
            if t['test_id'] == k or t['test_id'].startswith(k + ' '):
                t['before'] = short(t.get('final'))
                t['final'] = {'verdict': v, 'result': r, 'binary': EXBIN.get(k, BIN)}
                t['X'] = t['T'] = None
                t['notes'] = 'before: %s %s' % (t['before']['verdict'] or '', short(t['before']['result'], 110)['result'])
                if k in EXNOTE:
                    t['notes'] = EXNOTE[k] + ' ' + t['notes']
    R['extra_tests'] = [x for x in R['extra_tests'] if not x['test_id'].startswith('2-source')] + NEW_EX
    fill(R, sub)
    R['title'] = 'IDORT (Ma+2026, arXiv 2503.16627) + Jiang 2012/2021 test suite: FIN3 final scheme (fin3-1012 281dd686)'
    R['date'] = '2026-10-13'
    R['lede'] = ('Every IDORT (Ma et al. 2026, arXiv 2503.16627) test we run, plus the Jiang 2012/2021 tests where IDORT has none, '
                 'rerun with the FIN3 final AthenaK rad_m1 implicit VET scheme (fin3-1012 281dd686: restarts bitwise, '
                 'implicit_predictor_frozen corrected to false). Every figure is regenerated from the FIN3 runs; each test shows '
                 'the Cartesian (vet_sc) result and its spherical-wedge (vet_gd, production path) counterpart, with the previous '
                 'result as a short "before" line. No PENDING rows; FIN3 does not pass all tests: 8 OPEN: test 17 spherical side (rswv T), test 20 spherical side '
                 '(sp Marshak front), rw_p10t1 T ladder, T9v rswv T, T9tv rswtv, rlatt / rlatt_om5, T5 Marshak T ladder (marsh_r1_T), '
                 'single opaque disc. Final count: PASS 50 / OPEN 8 / PARKED 4 / EXCLUDED 2 / PENDING 0.')
    R['badge_note'] = ('Badges show the FIN3 final verdict (fin3-1012 281dd686); "Before" is the previous page state (FIN, '
                       'method-1012 78b2c961, with the d3bf42d7 jse fix and the twinfix candidate for T7). PENDING = final-config '
                       'rerun still running (Caltech / DeltaAI / Viper chains, Raven queue). Reproduced (qualitative) counts as PASS; '
                       'known-open items as OPEN; not-feasible tests as PARKED; EXCLUDED = kept outside the suite by user decision (multi-source bubbles, 2-source retardation), not counted, numbers for information only.')
    R['before_label'] = 'Before (FIN 78b2c961)'
    R['binary'] = {
        'final': 'fin3-1012 281dd686 = 1fdbe56e (FIN3 + bitwise restarts + output cadence) + implicit_predictor_frozen = false '
                 '(fresh default); binaries fin3_1012/bin/athena_{cpu,he_cpu,cpu2k}_1fdbe56e / _914e5433 with the key named',
        'sources': 'Cartesian: IN\'s final-config reruns fin3_1012/runs_281dd686 (cart281); sp: fin3_1012/sp_FIN3_281dd686 '
                   '(= the 914e5433 FIN3 arm + the 1fdbe56e two-stage reruns; implicit_predictor_frozen is bitwise-inert on sp); '
                   'P2: p2_FIN3_281dd686; T7 432 dirs: idort_page/reruns/fin3/sp_shd_FIN3_432',
        'previous': 'method-1012 78b2c961 (FIN) and bsgcand-1012 d3bf42d7 (jse fix); prodfinal-1011 4253af30 (PRODd) before that',
        'latest_candidate': BIN}
    R['scheme'] = {
        'name': 'FIN3 final = FIN2 (FIN + implicit_thin_js_face hcface, vet_gd_rotate_slope correct) + vet_gd_twin_state trans + '
                'vet_rotate_avg 4 + UN\'s fresh defaults (vet_interp auto, vet_col_lat_cut_thin 1, vet_bc_inner / vet_bc_reflect '
                'auto, vet_source_noesrc, implicit_predictor_frozen false) + BB\'s vet_sc rotation on the shared clock and '
                'transverse marshak BCs (implicit_bc_x2/x3) + bitwise restarts',
        'final_keys': [
            ['time_scheme', 'hesdirk2'], ['time2_vet_sc / time2_vet_col', 'predict'], ['time2_vstage', 'true'],
            ['implicit_vimp', 'true where gas is live'], ['time2_pred_cap', '0.5'], ['time2_start', 'be'],
            ['implicit_gas_src', 'conservative2'], ['implicit_cs2_hll', 'drag'], ['implicit_cs2_shift', 'false'],
            ['implicit_vimp_bc', 'true'], ['implicit_flux / implicit_flux_faces / implicit_flux_beam', 'blend / all / halfrange'],
            ['implicit_blend', 'idort_f (tau0 1, alpha 1, flo 0.3, fhi 0.6), xthin 0'],
            ['implicit_thin_switch / corr / hcsf', 'fse / hc / j'], ['implicit_thin_recon', 'jse'],
            ['implicit_thin_js_gd / implicit_thin_js_face', 'x1 / hcface (FIN2)'], ['implicit_thin_bq / bq_gd', 'lin / cell'],
            ['implicit_thin_hcface / hcface_gd', 'tau / x1'], ['vet_seg', 'jump'], ['vet_bc_diff', 'face'],
            ['vet_gd_symref', 'col (vet_gd)'], ['time2_vet_col_source', 'blend (vet_col / vet_gd)'],
            ['vet_gd_rotate_dt / align / slope / stretch', '~20 x typical dt / true / correct (FIN2) / 0.02'],
            ['vet_gd_twin_state', 'trans (FIN3, SH)'], ['vet_gd_twin_weight', 'none'],
            ['vet_rotate_avg', '4 (FIN3, AV; vet_gd and vet_sc)'], ['vet_sc rotation', 'on the shared clock (FIN3, BB)'],
            ['implicit_bc_x2 / implicit_bc_x3', 'marshak transverse BCs where named (FIN3, BB)'],
            ['vet_interp', 'auto (FIN3, UN #6)'], ['vet_col_lat_cut_thin', '1 (FIN3)'],
            ['vet_bc_inner / vet_bc_reflect', 'auto (FIN3, UN #2)'], ['vet_source_noesrc', 'fresh default (FIN3, UN #11)'],
            ['implicit_predictor_frozen', 'false (FIN3 final, 281dd686)'], ['implicit_hr_face_gd', 'none (deferred)'],
            ['restarts', 'restart = continuous, bitwise (FIN3 final)']]}
    R['final_pass'] = {'binary': BIN, 'arm': 'FIN3 (Cartesian runs FIN3_*, sp arm HFfZ + FIN3 keys)',
                       'note': '/raven/ptmp/jinma/fin3_1012/NOTE-fin3.md; /raven/ptmp/jinma/idort_page/NOTE-fin3page.md',
                       'date': '2026-10-13'}
    import collections
    cnt = {}
    for s in ('cartesian', 'spherical'):
        c = collections.Counter()
        for t in R['tests']:
            d = t.get(s)
            if isinstance(d, dict) and isinstance(d.get('final'), dict):
                c[d['final']['verdict'].split()[0]] += 1
        cnt['tests_' + s + '_FIN3'] = dict(c)
    cnt['extra_tests_FIN3'] = dict(collections.Counter(t['final']['verdict'].split()[0] for t in R['extra_tests']))
    cnt['spherical_extra_FIN3'] = dict(collections.Counter(t['final']['verdict'].split()[0] for t in R['spherical_extra']))
    R['counts'] = cnt
    R['verdict_rules'] = ('final = FIN3 final scheme (fin3-1012 281dd686) against the page criteria; before = the previous page '
                          'state; PENDING = final-config rerun still running (never old numbers as final).')
    ch = {'date': '2026-10-13', 'worker': 'PG (FIN3 page), from IN / SH / BB / AV / RT / UN',
          'note': '/raven/ptmp/jinma/idort_page/NOTE-fin3page.md; /raven/ptmp/jinma/fin3_1012/NOTE-fin3.md',
          'items': ['all tests and figures regenerated with the FIN3 final scheme (fin3-1012 281dd686: restarts bitwise, '
                    'implicit_predictor_frozen false); previous results kept as a short "before" line',
                    'T7 spherical shadow: twin_state trans now in the scheme; 432-direction FIN3 rerun',
                    'P2 in the new definition (wide channel pinc_ha20_sr30, AV\'s p2eval): PASS',
                    'bubbles: outflow BCs, marked retardation-limited (RT: the time-dependent reference reaches 0.213 / 0.214 / 0.160)',
                    'cold opaque spherical shell: fixed in FIN3 (front cell 1.000, was 28-211x)',
                    'user-approved TEST criterion change (OD): opaque-disc reference = the staircase exact (the discrete disc the grid holds), also for the bubbles; single opaque disc FIN3 0.784 / 0.873 / 0.927 at n256 / 512 / 1024 (first order), OPEN (M1 corner deficit); bubbles on the staircase reference L1 0.214 / 0.249 / 0.135 (ring 100; information only, test excluded); candidate vet_path = cells (0.942 at n512) shown, NOT adopted',
                    'user decision: the multi-source bubbles test (IDORT 4.2.2, test 12; ring 100 / 200, Gauss, periodic, bubu) and the 2-source retardation row (the deferred light-travel item) are kept OUTSIDE the test suite: new badge EXCLUDED, not counted as PASS / OPEN / PENDING; their numbers stay as information only',
                    'user-approved TEST changes (OT, open-1013): rswv deck change (implicit_f0_init = zero + outflow gas top, rswv decks only, code defaults unchanged) fixes rswv X (1.73-1.95, interior 2.0); rswv T stays OPEN (v / E / F1 1.2-1.3, same per-step term as rw_p10t1)',
                    'user-approved TEST criterion (OT): rlatt_om5_T judged on the finer pairs (levels 4-16): E and F1 converged there (E 2.04 / 2.46, F1 2.05 / 2.30), but the row stays OPEN (user decision): v_phi 1.61 / 1.49 on the finer pairs and dens 1.73 at the finest pair; OT is checking whether v_phi is the same per-step conservative2 velocity term. The rlatt X F1 floor and rswtv stay OPEN',
                    'IN final patch: test 5 (hcyl 512^2), 10 / 11 (dd), 18 (shocks M 2 / 3 at 1024, M 3 at 512), 21 (pipe) and G5 go from PENDING to PASS; no PENDING rows left',
                    'test 21 pipe caveat: PASS on the end state, but the probe-5 arrival time is 4.20 / 3.31 (c dt/dx 220 / 2200) vs FIN2 3.81 / 3.71, so FIN3 is less dt-consistent than FIN2 there; the criterion is our own',
                    'final count PASS 50 / OPEN 8 / PARKED 4 / EXCLUDED 2 / PENDING 0; FIN3 does not pass all tests (8 OPEN)',
                    'Cartesian short runs were made with the Viper 1fdbe56e CPU build (bitwise to 281dd686 there); a Raven 281dd686 rerun reproduces every page metric to print precision, and only Picard counts / round-off floors changed (rows 1, 2, 7)']}
    R['changes'] = [ch] + [c for c in R.get('changes', []) if c.get('date') != '2026-10-13' or 'PG' not in c.get('worker', '')]
    json.dump(R, open(P + '/idort_results.json', 'w'), indent=1)
    print('wrote idort_results.json', R['counts'])


if __name__ == '__main__':
    import sys
    sub = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    main(sub)
