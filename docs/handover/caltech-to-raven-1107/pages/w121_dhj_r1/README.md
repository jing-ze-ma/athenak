# Hot Jupiter Round 1 page (WASP-121b incl. r1m30) — https://claude.ai/artifact/A5omHrtMwF5roVtjNjDqq7

This is the only page that shows the WASP-121b round-1 runs; there is no separate r1m30 page. Last READY 11-06 17:38:
w121h 302.4, w76 301.6, toi1518 301.3, w178 159.1, w189 134.1, w121m100 93.0, w121m30 116.4, w121m7 16.6, w121m8 89.1 rot (of 300).
Counter ~v3+ (11-02); last publish 11-06.

## Files (unchanged copies)

- `refresh.sh [status]` (from `$P/dhj_r1/`): quick `python3 ana.py html` (2-core CPU step in a held alloc: candidates 5274196 /
  5719508, else login), then a DETACHED catch-up `ana.py full` (NPROC=4, `srun --overlap ... -c 4 --mem=40G --gres=none`) that writes
  `READY` (line 1 + files: index.html + img/*.png); also launches `fe_step.sh`.
- `fe_step.sh` — optional Fe I high-resolution blueshift update (viper's fe_hires pipeline in `$J/fewind_1102/fehires_fast`,
  venv_fehgpu, picks a free test GPU); the page builds without it. Static outputs copied in `static_inputs/fehires/`.
- `tools/ana.py` — reductions (cached per dump), figures, index.html. RUNS: w121h `w121r1_1101/runs/p1`, w76/toi1518/w178/w189
  `r1planets_1101/<planet>/runs/p1`, plus every `w121r1_1101/runs/r1*m*` dir with dumps and `run.ia*.log` that is not STOPPED
  (-> w121m30, w121m100, ...). Grey overlays from old r1b runs (OLDCACHE). EOS `w121r1_1101/anchor/eos_1x.txt`, NROT 300.
- `tools/budget.py` (conservation budgets of the hydro p1 runs), `rochedhj_ana/dhjcs.py` (cubed-sphere dump reader/geometry; imports
  bin_convert from `/resnick/home/jingze/ATHENAK/athenak/vis/python`).
- `READY` (current file list, 11-06).

## Inputs

hst `dhj.user.hst` / `dhj.hydro.hst` (or `dhj.mhd.hst`) of each run, `bin/*.bin` dumps (per-dump scalars + profiles on
log10 p = -6..2 bar, levels 1e-3 and 1e-1 bar), run logs `run.ia*.log`, decks. Static: `$J/obscheck_1102/page/` (phase curves,
emission, time-resolved Fe; copied to `static_inputs/obscheck_1102_page/`), fehires_vs_seidel.{png,txt}.

## Caches (committed here, small)

`cache/dhj_p1_page_cache/` (3.3 MB; per-run scal.json/profiles + budget/summary.json) and `cache/dhj_r1b_page_cache/` (1.4 MB, grey
overlays). Copy them to the paths in ana.py (`CACHE`, `OLDCACHE`). The r1m30 history is in `dhj_p1_page_cache/w121m30`.

## Publish

`index.html` + every `img/*.png` listed in READY (no inlining), url explicit.

## Machine-specific bits

`J`, `PAGE`, `CACHE`, `OLDCACHE`, `EOSF`, `R1`, run dirs, `OB`/`FE` static paths in ana.py; `sys.path` in ana.py and dhjcs.py; alloc
ids/nodes in refresh.sh; `module load python/3.11.6-gcc-11.3.1-rphh4kv py-matplotlib/3.8.2-gcc-11.3.1-7d3v7ui`.
