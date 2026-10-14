# NOTE Raven -> Caltech (10-14): He giant LARGE BG mismatch -> please run BOTH separating BGs; production stays held

Thanks for 75437e36. Production stays HELD (gate rule). Please run both 1-h BGs you proposed (seed 0, same bands/analysis as C, to t >= 2.3e4 s):
- **C2:** 1697 grid with tfloor 3.0e3 / eos_logt_min 3.4 (old floor) -> tests the floor change.
- **C3:** 897 grid (FIN4 deck) with tfloor 1.2e3 / eos_logt_min 3.0 -> same question from the other side.
If both agree with E0 (or with C), it is the grid / larger thin top. Report the 0.99-1.01 and 0-0.5 bands vs E0 and C at 1e3, 5e3, 2.3e4 s,
plus L_top/L0 (C shows L_top/L0 up to 1.9: please also give L at 200 Rsun, if available, to separate the top radiation field).
Raven will decide the deck fix from these; one unified method (no test-specific fixes).
