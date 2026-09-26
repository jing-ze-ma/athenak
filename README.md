# AthenaK run data, 2026-09-26 (viper -> Caltech handover)

Data branch only: it has no code and shares no history with rt-integration. Code: branch `rt-integration`;
start at `docs/handover/HANDOVER-2026-09-26.md` there.

`athenak_data_2026-09-26.tar.gz` (md5 057e168a9a1559e9eb419bf362dae0b1) unpacks to `athenak_data/`:
- `exo_fms_ck/`: correlated-k tables for 1x solar metallicity, placeholder `ATHENAK_CK_DATA`. It is also the repo's gitignored `data/exo_fms_ck`.
- `ckdata10/`: the same at 10x solar, placeholder `ATHENAK_CK_DATA10`.
- `he_box/`: He-box M1 IC, opacity and a_rad files, placeholder `HE_BOX_DATA`.
`MANIFEST.md5` lists the per-file checksums.

Download without cloning the code:
    git clone --depth 1 --branch data-2026-09-26 --single-branch https://github.com/jing-ze-ma/athenak.git athenak-data
or
    curl -LO https://github.com/jing-ze-ma/athenak/raw/data-2026-09-26/athenak_data_2026-09-26.tar.gz

Provenance and terms: the Exo-FMS tables come from https://github.com/ELeeAstro/Exo-FMS_column_ck (Lee et al. 2021,
MNRAS 506, 2695). The upstream repository has no licence file, so these copies are for research use only. Cite Lee et al. (2021)
and the underlying line lists, and check with the author before publishing. The hiT/hiT2 and 10x tables are derived
from them (data/exo_fms_ck/HITEMP.md, tools/gen_hitemp.py on rt-integration). See data/exo_fms_ck/PROVENANCE.md.
