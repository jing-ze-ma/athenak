# NOTE from Orion: low-density Planck/Rosseland task DONE (re TASK-2026-10-09-orion-lowrho-planck)

Results: branch `agcar-opac-orion-1009`, commit b9c9e4f6. Read `docs/handover/NOTE-2026-10-09-orion-lowrho-planck.md`.
Tables are on the ext2 nodes (log T 3.45-4.3 step 0.025, log rho -21..-12 step 0.05), X 0.36 Z 0.02 GS98:
`docs/handover/agcar-opac-orion-1009/tables/{planck,rosseland}_lowrho_orion_gs98_x0.36_z0.02.txt`.
These are LTE means built from the Kurucz line lists, with FastChem chemistry and Doppler line profiles.

## Validation (orion minus reference, dex)
- Rosseland vs Ferguson 2005: median 0.000, max 0.08. vs TOPS: median -0.013, max 0.10.
- Planck vs Ferguson 2005: median +0.025, p10/p90 -0.37/+0.30. vs TOPS: median -0.28, p10/p90 -0.53/+0.06.
- Planck by temperature:
  - Good at T <= 3500 K and 6300-7900 K.
  - At 9-10 kK, Ferguson and TOPS disagree with each other by up to 1 dex; orion lies between them.
  - At 11-20 kK, orion is 0.3-0.55 dex LOW against both references. This is unexplained: it is not the predicted
    Fe-group lines, the hydrogen n_max, the continuum, or the newer C III line list. Treat orion's kP above 11 kK as
    a lower bound, short by up to about 0.5 dex.

## ext2 extension error at rho 1e-20..1e-15 (log T >= 3.5)
- Planck: errors from -1.0 to +0.8 dex. kP is non-monotonic in rho: it drops about 1 dex at each LTE ionisation of
  Mg, Fe, Ca, Si or C, so no power law can follow it.
- Rosseland: the ext2 extension is fine (median 0.00, max 0.12 dex). Keep it.

## Proposed stitch (Orion's proposal; the user decides)
- Planck, below each source floor rho_f(T) that ext2's power-law extension starts from:
      kP_new(T, rho) = kP_ext2(T, rho_f(T)) * kP_orion(T, rho) / kP_orion(T, rho_f(T))
  This takes the density shape from orion and the normalisation from the source table ext2 already uses at that T.
  It is continuous at the floor by construction, leaves every real Ferguson/TOPS cell unchanged, and absorbs the
  11-20 kK normalisation deficit.
- Outside orion's T range (log T < 3.45 or > 4.3): keep ext2, and ramp the correction factor to 1 over one or two
  T nodes at each edge.
- Rosseland: keep ext2 unchanged.
