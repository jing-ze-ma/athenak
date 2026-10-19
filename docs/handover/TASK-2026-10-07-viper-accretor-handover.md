# TASK for viper (from DeltaAI, user 10-07): hand the accretor sim over to DeltaAI + one physics question

**User 10-07:** the accretor (ry_per_accretor / Plaskett) simulations now run on DeltaAI and all decisions about them
are made there. Please hand over everything you have and answer the question in section 2. Reply as
`docs/handover/NOTE-2026-10-0x-viper-accretor-handover.md` on this branch (`accretor-1006`, fast-forward, never force).

## 1. Hand-over: please push to `docs/handover/accretor-handover/` on accretor-1006

Everything you thought about this problem, so DeltaAI can continue without asking again:
- **Design docs** verbatim: `/viper/ptmp2/jinma/accretor_1006/DESIGN.md`, `PLASKETT.md`, and any other design / trade-off
  / results notes for ry_per_accretor (the inputs cite `# Design / trade-off: /viper/...`).
- **Parameter derivations**, with numbers: m_acc, m_don, a, P, r_acc, t_don, mu_don, env_npoly, **env_rho_ph = 0.3**,
  env_cs_ph, env_cs_stream, env_cs_amb, env_amb_rho, env_amb_k, t_relax values, sponge keys, spin 1.0 vs 7.2, env_r_spin,
  r_in / r_out (0.85 d_L1), the grid stretch, wb_rmax, the stream window (nsig 3) and the stream_width default.
- **Run history**: every arm so far (q10_s1 dt collapse at 0.077 orbit and its diagnosis, q12_a1 / q12f_a1, q12_s45, smokes,
  anything else), with what each one taught you and why the keys changed (6588011b, ea04cd19, ...).
- **Open problems and planned next steps**, including what you would run next and why; known pgen limitations
  (e.g. thermo = adiabatic only with the envelope mode; isothermal-only paths).
- **Analysis scripts** beyond page/mkfigs_plaskett.py, and the scientific goal of the project (which quantities matter:
  j of accreted gas, spin-up, Mdot, ...), plus any literature you relied on.
- Paths on viper of the run directories (DeltaAI cannot pull from viper; small files in the repo please, big ones list only).
- Your pending viper jobs 12116831 (q12_s1) / 12116832 (q12f_s1): DeltaAI's copies of both arms ran to (near) tlim today
  (NOTE-2026-10-07-deltaai-plaskett.md); state what you did with yours.

## 2. Question: how was env_rho_ph = 0.3 (and the stream density scale) derived?

DeltaAI plans to set the stream width from Ryu et al. 2025 (A&A, arXiv:2505.18255, isothermal L_in overflow, which matches
our relaxed-to-c_ph^2 stream). For this binary (q = 1.125, c_s,don 20.96 km/s, Omega 1.9708e-5 /s):
- L1 coefficients A = -16.99, B = 6.99, C = 7.99 (exact; the paper's A(q) fit gives -16.79); d_L1 15.884 Rsun = pgen.
- Analytic isothermal stream at L1: sigma_y = c_s/(Omega sqrt B) = 0.578 Rsun, sigma_z = c_s/(Omega sqrt C) = 0.541 Rsun.
- Their q = 1 isothermal SIMULATION profile (digitised from their Fig. comparison): FWHM 0.932 x analytic, peak 0.935 x,
  peak shifted 0.13 c_s/Omega to the trailing side; their fits vary only ~1e-3 between q = 1 and 1.125.
- => proposed `problem/stream_width = 0.54 Rsun` (now default c_s/Omega = 1.529 Rsun, 2.8x wider; measured FWHM in
  q12_s1 at r_out 3.5 Rsun as imposed). Tilt: paper -25 deg vs our ballistic launch ~28 deg at r_out: consistent.

Narrowing the stream at fixed Mdot raises its peak density by ~ (1.529^2)/(0.54 x 0.50) ~ 8.6. **Question:** was
env_rho_ph = 0.3 ("rho_ph / stream peak at Mdot 1e-4", PLASKETT.md) computed from Mdot = 1e-4 Msun/yr with a stream
cross-section of 2 pi (c_s/Omega)^2 and which stream speed? Please give the formula and numbers, and say which of these is
the physical input you intended to hold fixed: (a) Mdot = 1e-4 (then env_rho_ph -> ~0.035, and env_amb_rho / sponge keys
rescaled?) or (b) the density contrast 0.3 (then Mdot ~1.2e-5). Anything else in the setup that assumed the old width?
