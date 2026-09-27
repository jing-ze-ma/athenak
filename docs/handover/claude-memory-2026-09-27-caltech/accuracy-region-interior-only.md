---
name: accuracy-region-interior-only
description: Massive-star boxes/wedges - the atmosphere near the top boundary does not matter (like p < 1e-6 bar for dhj); judge accuracy and cost below it
metadata:
  type: feedback
---
User, 09-26: "the atmosphere especially near the top boundary doesn't matter for massive stars ... similar to above 1e-6 bar for dhj."
**Why:** the science is the interior and envelope (convection, the Fe opacity bump, transport). The thin top is a boundary buffer.
**How to apply:**
- Accuracy gates, tolerance sweeps and A/B tests on He/massive-star runs: measure deviations below the photosphere only (Rosseland tau >~ 1, excluding the top few cells and the sponge). Use the flux through the interior or photosphere, not F1top at the boundary.
- A cheap or crude treatment of the top (closure, tolerance, resolution) is acceptable if the interior is unchanged.
- CAVEAT (user, same day): the top must still not significantly limit the timestep, crash the code, or use much computation. It is irrelevant for accuracy but NOT for cost or stability. Report the dt limiter location, failures and iteration counts in the top separately.
- The dhj equivalent is [[resolution-only-below-1e-6-bar]]. See also [[no-accuracy-sacrifice]].
