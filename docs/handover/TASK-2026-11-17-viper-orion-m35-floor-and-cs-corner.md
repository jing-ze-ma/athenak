# TASK 2026-11-17 viper -> Orion (USER): (1) M35MMS floor runaway: would stag2 avoid it? (2) cubed-sphere corner unit test for RSG

## 1. M35MMS: the conservative gas write-back pumps floor gas (production rad_m1)
- **Run:** Will Schultz's 35 Msun MS envelope. 3-D spherical-polar wedge, r 9.7-15.3 Rsun, general EOS, MESA Rosseland + TOPS Planck, production
  binary 0940c7d1_t9fix2 (prod-1118 line), `rad_m1/implicit_gas_src = conservative2`, hesdirk2 + time2.
- **Symptom:**
  - From t = 1.645 d the floor gas above the photosphere (rho 1e-15 .. 1e-19, r > 14 Rsun) heats to the EOS ceiling of 1e8 K.
  - v_r there reaches the 1e9 cm/s vceil; dt falls 7 -> 0.35 s.
  - The first collapse cell is the domain corner (first active theta cell, second phi cell, r = 14.3 Rsun).
- **Trigger in the log:**
  - At cycle 18159, time2_inc_kecap (INCKECAP) fires on 3 floor cells with KE_not_given 9.4e31 erg, i.e. a stage momentum increment dm/rho ~ 1e12
    cm/s.
  - From cycle 18264 on, BE solves "kick cells to |v| > 1e9 (max 2.79e12)".
  - A 2-rank rerun from 1.5 d does not reproduce the first event (trajectory-dependent).
- **Mechanism (pump proven, trigger hypothesis):**
  - conservative2 gives the gas minus the radiation's own change, dm = -(F'_c - F^old_c)/(c chat) - sum Z_f/2, and likewise for energy.
  - In floor gas F/(c chat) ~ 1e-7 against rho 1e-15 .. 1e-19, and E_rad ~ 1e4 against e_gas ~ 0.02 erg/cm^3.
  - So any small residual or mismatch of the solve becomes dv ~ 1e6 .. 1e12 and a x500 jump in e/rho.
- **Tests (2 ranks, MI300A):**
  - Recovery from the post-event rst: conservative2 keeps 427 vkicks and dt 0.39 s. implicit_gas_src = matched (with implicit_lzx = false, which
    matched requires) gives 1 vkick and dt rising 0.49 -> 0.78 s.
  - Forward from a healthy 1.644 d rst through 1.70 d with matched: clean, dt 7-8 s, energy and L_top equal to conservative2.
  - Cost: L_z is no longer exact (not measured).
- **Production decision (USER, 11-16 ~14:50 PST):** M35MMS is restarting from 1.5 d with matched + lzx off. The proper fix proposed: matched (or a per-cell
  exact local exchange, as clean_yres = local / clean_gas_exact do) only in floor cells, with a physical switch on e_gas/E_rad (~2e-6 in the floor,
  > 1e-2 in the star).
- Branch: fork m35-floorfix-1116 (diagnostic key dbg_inckc_ke; optional floor T ceiling `<hydro>/tceil_kelvin`).
  Write-up: viper /viper/ptmp2/jinma/m35floor_1116.

**USER's question to you: can this be avoided with the new staggered scheme (stag2)?** Please answer from stag2's design and, if cheap, a test:
- How does stag2's gas update take the radiation's momentum and energy? (Your NOTE said "the gas gets only physical quantities: the force h chi (F - A)/c
  and the exchange on the solution; the linear-solve residual stays in E".)
- Is the floor-gas amplification (residual / tiny rho) therefore absent by construction? Are there other routes: the F solve's residual, the energy
  write-back, positivity fixes, KE caps?
- A cheap CPU check: a 1-D or 2-D sp atmosphere with a floor of 1e-15 .. 1e-19 above a photosphere and E_rad/e_gas ~ 1e6 in the floor, run with the
  cell layout (clean_gas_exact on/off) and with stag2. Look at max |v| and T in the floor over a few hundred cycles. Compare with production rad_m1
  conservative2 vs matched if you can.

## 2. RSG E cube-corner jet: hydro-only cubed-sphere corner unit test (USER: go)
- RSG_BLOB_1116 (viper /viper/ptmp2/jinma/rc18_1110/RSG_BLOB_1116.md): the 0/1/5 corner jet came back in the rerun from 6988 d, within 10-25 % of the
  old branch. The keys are cs_seam_rho_guard = 4, cs_seam_blockend = extend, cs_vertex_3cell = true, cs_vertex_mono = true, cs_cross_mono = true.
  Caltech used guard = 0 and ran 1822 d with no corner event.
- **Please run on CPU:**
  - hydro-only cubed sphere;
  - a smooth overdense blob (or an outward velocity pulse) placed (a) at a cube corner, (b) at a panel centre, (c) on a panel seam away from corners;
  - arms: production keys as above; cs_seam_rho_guard = 0; cs_vertex_3cell = false; both off.
  - Measure the anisotropy of the outgoing front: corner vs centre speed and amplitude, and any focusing or jet at the corner.
- Reply as NOTE-2026-11-17-orion-viper-cs-corner.md. viper will run GPU twins of RSG only if the CPU test points at a key.
