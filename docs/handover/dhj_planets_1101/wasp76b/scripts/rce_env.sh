# environment for rce_w121.py / rce_f5.py / transit (wasp76b, 1x solar: ck/CK1X tables + /viper/ptmp2/jinma/deepconv_0925/dump EOS dump); source, then set RCE_AP, RCE_G0, RCE_OUT
export RCE_TEQ=2215.6 RCE_TSTAR=0 RCE_MET=1x RCE_ALB=0.20 RCE_RESTOL=1e-7
export RCE_SWFILE=/viper/ptmp2/jinma/dhj_planets_1101/wasp76b/ck/CK1X/sw_flux/sw_band_flux_W76_11.txt
export RCE_CKDIR=/viper/ptmp2/jinma/dhj_planets_1101/wasp76b/ck/CK1X
export RCE_EOSDIR=/viper/ptmp2/jinma/deepconv_0925/dump
export TR_TEFF=6329
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
