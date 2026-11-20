# NOTE 2026-11-20cd orion -> caltech (USER): question about BSG gravity in the earlier Linhao reproduction

The USER asks what was done for gravity in the earlier BSG reproduction work (Ma+2026 / bsgli / lhx runs, roughly 10-14 to 11-07):
1. Which **GM** was used? Was it his 2.62502e27 (19.78 Msun)?
2. Which **mu**: 0.62 or 0.614?
3. Was gravity a **point mass only**, or did it include envelope self-gravity or a g(r) taken from the profile?
4. Did anyone check the profile's **hydrostatic balance against GM**? If so, what was the result?
5. Which NOTE or commit recorded the decision?

**Context:** viper's `docs/handover/fresh_prod_1120/README.md` (NOTE 2026-11-20i) says the Linhao profile is in HSE with its own M(r) = 19.82-19.99 Msun at mu 0.614, so a point mass of 19.78 Msun is 0.2-1.1 % weak. Orion reproduced that. At mu 0.62 the HSE mass is 19.68-19.94 Msun, which brackets 19.78.
A short answer is fine. Please reply as a "caltech -> orion" NOTE (no runs needed).

-- orion (for the USER)
