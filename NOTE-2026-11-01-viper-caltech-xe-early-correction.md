# NOTE 2026-11-01 viper -> Caltech (USER): first in-code x_e tracer result, and a CORRECTION: the tracer does NOT remove the binding resistive dt

Our 2-rotation x_e tracer run on the s10 (10 G, rot 320) restart has reached 0.8 rot. Analysis: viper /viper/ptmp2/jinma/xetracer_1101/analysis/F_early.md. The binary is the z-form tracer (61511b84), without rkg_lean.

## 1. Night electron fraction: the tracer works
Median F = x_e,tracer / x_e,Saha on the night side (zenith > 120 deg), against the offline estimate:

| log p [bar] | -6..-5 | -5..-4 | -4..-3 | -3..-2 | -2..-1 |
|---|---|---|---|---|---|
| tracer (0.8 rot) | 1400 | 4200 | 590 | 51 | 1.2 |
| offline | 12 (248 with vertical motion) | 315 (523) | 780 | 83 | 1.2 |

- The night eta_O at 1e-5..1e-3 bar drops by 3-4 orders of magnitude (1.5e14 -> 5e10 cm^2/s at -5..-4).
- The night Rm goes from 0.09-0.35 to 28-1100.
- The terminators and the day side stay at F ~1 below 1e-6 bar, as expected.

## 2. CORRECTION to our earlier NOTEs (651f6b733, 5e4f80389)
We wrote that with the tracer "the stiff layer mostly goes away". **It does not.**
- The global binding h^2/eta sits at 1e-2..1e-1 bar (0.249 s Saha -> 0.253 s tracer). There x_e stays Saha, because recombination is fast at that density.
- RKG stage counts: 32-34 per cycle with the tracer vs 34 without, i.e. about 3-6 % fewer stages.
- So the r1m RKG cost is essentially unchanged.
- The tracer is a physics fix: night induction and magnetic drag at the Fe level. It is not a cost fix. Please keep that in mind for the implicit-resistivity scoping: the binding layer is deeper than our earlier estimate implied.

## 3. A defect found (being fixed): floor dilution at the night top
- At 1e-8..1e-6 bar on the night side, F = 0.01-0.1: x_e falls below Saha and eta_O rises 4-18x.
- **Cause:** the Alfven floor adds mass with zero electrons at every conversion.
- **Fix:** the floor mass inherits the cell's electron fraction. It goes on fork/tracers-defaults-1101 next. It matters only in the floored top.
