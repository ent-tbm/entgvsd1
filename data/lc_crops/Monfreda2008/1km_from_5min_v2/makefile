# This makefile processes the original 5-minute Monfreda .asc files,
# found in ../5min_orig/, with the following steps:
#
#  1. ../../entdata fetch:
#     Fetch ../5min_origdata/*.asc files from the NCCS Data Portal
#
#  2. asc_to_bin.f:
#     Sum up certain variables, store in GISS .bin format.
#     This is the original hntr4_monfreda2008_HXH.f, modified
#     slightly to:
#       a) Don't regrid
#       b) Variable names added to descriptors in .bin files.
#          (for giss2nc to pick up, see below)
#
# 3. giss2nc.cpp:
#    Convert output of (1) to NetCDF format.
#    (see https://github.com/citibeth/icebin/blob/develop/modele/giss2nc.cpp):
#
# 4. regrid_monfreda.F90:
#    Regrid output of (2) to 1km

all : Monfreda_crops_1km.nc Monfreda_crops_1km_norm.nc

ASC_FILES=../5min_origdata/C3C4/L1C3C4.asc ../5min_origdata/C3C4/L2C3C4.asc ../5min_origdata/form/L2form.asc ../5min_origdata/form/L3form.asc ../5min_origdata/form/L4form.asc ../5min_origdata/form/L1form.asc ../5min_origdata/legume/L1legume.asc ../5min_origdata/legume/L2legume.asc

# Fetch
$(ASC_FILES) :
	cd ../5min_origdata; ../../entdata fetch; gunzip `find . -name '*.asc.gz'`


Monfreda_crops_5min.bin Monfreda_crops_5min_norm.bin : $(ASC_FILES)
	gfortran -cpp -fconvert=big-endian -O0 -fno-range-check asc_to_bin.f; ./a.out

Monfreda_crops_5min.nc :
	giss2nc --input-file Monfreda_crops_5min.bin --output-file Monfreda_crops_5min.nc --endian=big

Monfreda_crops_5min_norm.nc :
	giss2nc --input-file Monfreda_crops_5min_norm.bin --output-file Monfreda_crops_5min_norm.nc --endian=big


Monfreda_crops_1km.nc Monfreda_crops_1km_norm.nc : Monfreda_crops_5min.nc Monfreda_crops_5min_norm.nc
	xent regrid_monfreda.F90; python3 copy_metadata.py

# -------------------------------------------------

deploy : Monfreda_crops_1km.nc Monfreda_crops_1km_norm.nc
	rsync Monfreda_crops_1km.nc Monfreda_crops_1km_norm.nc *_forplot.nc discover:/discover/nobackup/projects/giss_ana/pub/Ent_TBM/EntGVSD/data/lc_crops/Monfreda2008/1km_from_5min_v2/

clean :
	rm *.bin *.nc a.out fort.*

