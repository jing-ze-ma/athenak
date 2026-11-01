# environment for rce_w121.py / rce_f5.py / transit (WASP-189b); source, then set RCE_AP, RCE_G0, RCE_OUT
# W189_MET = 3x (default; CK3X + the 3x EOS dump) or 1x (repo data/exo_fms_ck + the 1x EOS dump)
export RCE_TEQ=2637.6 RCE_TSTAR=0 RCE_ALB=0.10 RCE_RESTOL=1e-7
export RCE_SWFILE=/viper/ptmp2/jinma/dhj_planets_1101/wasp189b/ck/CK3X/sw_flux/sw_band_flux_W189_11.txt
if [ "${W189_MET:-3x}" = 1x ]; then
  export RCE_MET=1x RCE_CKDIR=/viper/u2/jinma/ATHENAK/athenak/data/exo_fms_ck
  export RCE_EOSDIR=/viper/ptmp2/jinma/deepconv_0925/dump
else
  export RCE_MET=3x RCE_CKDIR=/viper/ptmp2/jinma/dhj_planets_1101/wasp189b/ck/CK3X
  export RCE_EOSDIR=/viper/ptmp2/jinma/w121prod_0929/met3/dump
fi
export TR_TEFF=8000
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
