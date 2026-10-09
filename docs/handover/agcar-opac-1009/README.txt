Copy of /viper/ptmp2/jinma/lbv_1008/agcar/tables_ext (10-09) for Orion (TASK-2026-10-09-orion-lowrho-planck).
The scripts keep their viper absolute paths. Map:
  /viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/scripts/          -> scripts/
  /viper/u2/jinma/ATHENAK/athenak/data/stellar_opac/planck_tools/ -> scripts/planck_tools/ (compare.read_ferg/read_repo, convert_tops.bilin)
  /viper/ptmp2/jinma/lbv_1008/agcar/tables/ (+ raw/tops_gs98_x0.36_z0.02.dat) -> sources/tops_only/ (TOPS rho 1e-14..1; bitwise region of ext2)
  /viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/raw/tops_lowrho/tops_lo.dat -> sources/tops_lowrho/ (TOPS rho 1e-21..1e-14)
  /viper/u2/jinma/ATHENAK/bench/m1_opac/ferguson05/{ross/g98.35.02.tron,ross/g98.5.02.tron,g98.pl.35.02.tpon,g98.pl.5.02.tpon} -> sources/ferguson05/
    (members of f05.gs98.tar.gz md5 96a1b1db / f05.g98.pl.tar.gz md5 69d8b6d3)
ext2 = scripts/build_ext.py fergR <outdir>. Not included: AESOPUS raw (only for the 'ext' variant), run1d/ outputs, png checks.
MD5SUMS covers every file here; MD5SUMS.tables_ext is the original viper MD5SUMS.
