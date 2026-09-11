---
name: validate-the-instrument
description: A gate that reports "nothing" usually means "I did not look" -- the recurring false-pass shapes on this project, with the specific Viper ones
metadata:
  type: feedback
---

**The single most expensive recurring mistake on this project.** A measurement that
reports "clean", "0", or "no output" is not evidence until the instrument has been shown
capable of reporting a failure. Seven distinct instances so far; five in one session.

**The catalogue of false passes seen here:**
- **Rank-local norms** faked an MPI bug (2943fd02), and a PRE-EXCHANGE array
  (`efld_resist`) faked an O(1) resistivity bug. Ask WHICH ARRAY a gate reads and AT
  WHAT POINT in the task list.
- **`mpirun: command not found`** -- Viper's MPI modules are hierarchical; without
  `module load gcc/14 openmpi/5.0` the whole battery is vacuous and prints "guards ran"
  (they refuse) and "1 of 1 dumps differ" (an unmatched glob).
- **Comparing a glob that matches nothing** -- "0 of 0 dumps differ" reads as a pass.
  `cubed_sphere_resist_smr.athinput` has **NO `<output>` block**, so it writes NOTHING.
  Add one and CONFIRM dumps exist before trusting any bitwise comparison.
- **`/tmp` is NODE-LOCAL on Viper.** A batch job whose `-o` path and inputs live under
  the session scratchpad fails in ~18 s on the compute node with no output file. Stage
  on `/viper/u2` (e.g. `bench/`).
- **`grep` silently suppresses matches in a file containing a NUL byte** ("binary file
  matches"), so a passing diagnostic reads as absent. Use `grep -a`.
- **`comm` on differently-sorted inputs** prints a wrong (empty) intersection while
  warning `input is not in sorted order` -- easy to miss, reads as "no violations".
- **`make | grep | head` kills make via SIGPIPE**, leaving a stale binary that then
  reports pre-fix numbers. And `git stash` from the wrong directory silently does
  nothing, making a before/after comparison vacuous.
- **A stale recorded baseline** is not a regression. When numbers move slightly, build
  the actual HEAD in a `git worktree` (safe -- no stash) and compare; the narrow-halo
  values were byte-identical and the memory entry was simply old.

**How to apply:** before believing a gate, make it FAIL ON PURPOSE, or check a negative
control. Prefer comparing DUMPS or an exact-solution norm over per-block intermediates.
When several alarms fire at once, suspect the harness before the code. See
[[measure-impact-before-claiming]], [[cubed-sphere-smr]].

## The working mpirun on a viper LOGIN node (2026-09-01)

A plain `mpirun -np 2` dies in `MPI_Init` with `Returned "Unreachable" (-12)`.  What runs:

```
export PATH=/mpcdf/soft/RHEL_9/packages/znver4/openmpi/gcc_14-14.1.0/5.0.10/bin:$PATH
OMPI_MCA_btl=self,vader,tcp OMPI_MCA_pml=ob1 mpirun -np 2 ./build/src/athena ...
```

The same prefix is what a fresh MPI build needs for `-D CMAKE_CXX_COMPILER=.../mpicxx`;
a bare `cmake -D Athena_ENABLE_MPI=ON` fails with "Package 'mpi-cxx' ... not found".
And a FAILED cmake configure still writes a cache: the next `cmake -B <same dir> -D
PROBLEM=...` silently kept `PROBLEM=built_in_pgens` and the binary refused the input.
`rm -rf` the build dir after any failed configure.

**2026-09-01, three more.**
* **`sacct` said COMPLETED, exit 0 -- in 38 SECONDS.**  The run had gone NaN; a NaN dt
  raced `time` to `tlim` and the code exited normally.  A clean exit status is not a
  clean run.  (Earlier the same trap wore a different mask: a cycle counter advancing
  while the history file already held a NaN row.)
* **`athena -d <dir>` APPENDS to an existing `.hst`.**  I re-ran an A/B into directories
  left over from a previous build and compared two binaries' output mixed together;
  it reported a spurious "DIFFERS" twice before I checked the row counts, which had
  GROWN (64 -> 128 -> 192).  Growing row counts in a fixed-length run are the tell.
  Clear the output directories between A/B rounds.
* **Overriding `output*/dt` on a RESTART silently does nothing** -- the output schedule
  comes from the restart file, so I got 3 dumps where I asked for 60 and briefly read
  that as "the failure is not being captured".  The way to get dense output near a
  failure is to **bisect `tlim`**: a dump is always written at tlim.

## Three process traps (2026-09-02 session), same family

* **`pkill -f <pat>` / `pgrep -f <pat>` match the Claude Bash wrapper's own command
  line**, so killing "the run" killed my own shell (exit 144) TWICE.  Select with
  `ps -o comm` or an explicit PID.
* **A waiter like `grep -q X <(tail -3 log)` false-positived** because the tail window
  still held the PREVIOUS resolution's line.  Count occurrences
  (`[ $(grep -c X f) -ge 2 ]`) instead of tailing.
* **A cs run's log looks stuck at cycle=0 and is usually fine** -- these runs print only
  at cycle 0 and at the end.  Run refinement arms in PARALLEL (256 cores, load ~5):
  ~15 min instead of ~2 h sequential.


## 2026-09-04: the BINARY-VINTAGE trap -- a gate on a binary that predates the feature

Four restart-based gates for `<mesh>/use_polar_average_b` came back "IDENTICAL" and were
each explained away (empty output dir, NaN vs NaN, zero steps, wall-clock limit). The FIFTH
explanation was the real one: `build_cpu_dhj` had last been built BEFORE the flag was added
(for the outer-x1 BC fix), and only the HIP and the two MPI test builds were rebuilt after.
The runs parsed the flag (`-n` dumps the MERGED INPUT, not what the code reads) and wrote
`use_polar_average_b = true` into their own restart headers -- while executing a binary
with no such code path. A one-step run on a freshly built binary changed the sample
pole-face cell 0.3885 -> 0.2496; the eight-step "gate" on the stale binary left it at
0.3865 with the flag on AND off.

**Rule**: before gating a new flag, prove the binary HAS it -- and NOT with `strings` or
`grep -a`: on this build they find NEITHER the new flag NOR `use_polar_boundary` in
mesh.cpp.o, a literal that is certainly there, so they gave a false "missing" twice and
aborted a valid gate. What works: (1) compare the object/binary mtime with the commit time
(`ls -la --time-style=long-iso ... mesh.cpp.o` vs `git log -1 --date=iso`); (2) a
FUNCTIONAL check -- one step with the flag on must change the state it targets (here the
sample pole-face cell 0.3885 -> 0.2496); (3) best, have the code print the flag from the
object that reads it. A parameter dump proves nothing about the executable.
