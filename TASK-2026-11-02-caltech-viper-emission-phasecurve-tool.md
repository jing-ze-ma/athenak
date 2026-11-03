# TASK 2026-11-02 Caltech -> viper (USER): 3-D emission / phase-curve tool for the dhj outputs?

**User request (11-02):** check the two finished 300-rot HYDRO dhj runs (WASP-121b, WASP-76b; cubed sphere C32,
1x, Roche tide, no drag) against observations other than Seidel 2026: JWST NIRSpec G395H / NIRISS phase curves,
Spitzer 3.6/4.5 um, HST WFC3 (WASP-121b); Spitzer 3.6/4.5 (May & Stevenson 2021) and JWST (WASP-76b).
Proper post-processing: 3-D ray tracing through the model with correlated-k opacities (not 1-D columns).

Caltech has no emission / phase-curve tool: no gCMCRT, no petitRADTRANS setup for dhj outputs, nothing in any
fork branch (searched all remotes for phasecurve/emission/eclipse/gcmcrt/petitrad).

## Please (reply in a NOTE on this branch or bsg-files-1009; Caltech does not block on it)
1. Do you have a 3-D emission ray tracer or a pRT/gCMCRT setup for AthenaK cs dumps (viper path, branch, inputs)?
2. If yes: push it (like fe_hires) or say how its opacity/chemistry is set up, so we can cross-check our numbers.

## What Caltech builds meanwhile (obscheck_1102)
- Ray tracer: parallel rays over the projected disc at each orbital phase through the 3-D cs grid (C32 x nx1 112),
  formal solution of the emergent intensity per wavelength x g-point.
- Opacity: petitRADTRANS 3 R=1000 correlated-k line opacities (H2O, CO, CO2, CH4, OH, HCN, TiO, VO, FeH, Na, K, SiO),
  premixed with random overlap on a (T, p) grid, abundances from FastChem (pyfastchem, 1x solar, gas phase + equilibrium
  condensation removal), H- bf/ff, H2-H2/H2-He CIA, Rayleigh; no clouds. Stellar spectra PHOENIX.
- Cross-check: the same tracer on the run's own 11-band premixed ck table (pkg/ck/CK1X).

Runs (Caltech, read-only): /resnick/groups/carnegie_poc/jingze/w121r1_1101/runs/p1,
/resnick/groups/carnegie_poc/jingze/r1planets_1101/wasp76b/runs/p1 (dumps rot 298-300).
