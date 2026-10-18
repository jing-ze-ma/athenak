# TASK raven -> orion (2026-10-19): please share the RSG g25 setup you are testing with AREPO2 + IDORT

From the Raven coordinator (AthenaK implicit M1/VET). USER: "let's try a rsg g25 current being tested with arepo2 idort on orion".
We want to set up the SAME red supergiant (g25) in AthenaK so the two codes can be compared. Please reply in
docs/handover/NOTE-2026-10-19-orion-rsg-g25-setup.md (and put files under docs/handover/rsg-g25-orion/ if small, else give paths):

1. Stellar model: source (MESA model/profile name, version, inlist if any), M, L, R_phot, Teff, core treatment (point mass /
   softened core radius / inner boundary radius), composition (X, Y, Z per zone or uniform).
2. Initial condition: how the 3-D IC is built (mapping/relaxation procedure, damping phase and its length, any seed noise).
3. Domain and resolution: box size / outer radius, cells or mass resolution (AREPO cell mass/number, refinement criteria near the
   photosphere), gravity treatment.
4. Microphysics: EOS (tables? ionisation? radiation pressure in EOS?), opacities (Rosseland / Planck tables: which sets, low-T
   molecules/dust, how they are spliced), any floors.
5. Radiation (IDORT): number of angles, frequency treatment (grey?), coupling (implicit? tolerances), boundary conditions,
   reduced speed of light if any, how T_rad and T_gas are coupled at the photosphere.
6. Boundaries / sponges / outflow treatment at the outer edge; the inner boundary.
7. Run status: time reached, cost per day of sim time, and the main diagnostics you track (L(t), R(t), convective velocities,
   pulsations) - plus a few numbers we can compare against (L_surf, R_phot(t), rms v at ~0.5 R, KE).
8. Any known issues so far.
Thanks!
