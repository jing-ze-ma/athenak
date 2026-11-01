# environment for rce_w121.py / rce_f5.py / transit (wasp178b, 3x solar: ck/CK3X tables + /viper/ptmp2/jinma/w121prod_0929/met3/dump EOS dump); source, then set RCE_AP, RCE_G0, RCE_OUT
export RCE_TEQ=2452.3 RCE_TSTAR=0 RCE_MET=3x RCE_ALB=0.15 RCE_RESTOL=1e-7
export RCE_SWFILE=/viper/ptmp2/jinma/dhj_planets_1101/wasp178b/ck/CK3X/sw_flux/sw_band_flux_W178_11.txt
export RCE_CKDIR=/viper/ptmp2/jinma/dhj_planets_1101/wasp178b/ck/CK3X
export RCE_EOSDIR=/viper/ptmp2/jinma/w121prod_0929/met3/dump
export TR_TEFF=9200
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
