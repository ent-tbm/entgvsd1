#!/bin/bash
#
set -e

for stage in \
    B01_bnu_laimax.F90 \
    B02_lc_laimax_modis_entpftrevcrop.F90 \
    B03_regrid_snowice.F90 \
    B04_veg_height.F90 \
    B05_carrer_mean.F90 \
    B06_albmodis_gridfill.F90 \
    B07_soil_albedo.F90 \
    B08_lc_laimax.F90 \
    B09_lc_lai_doy.F90 \
    B10_lc_lai_monthly.F90 \
    B11_reclass_annual.F90 \
    B12_reclass_doy.F90 \
    B13_reclass_monthly.F90 \
    B14_regrid.F90 \
    B15_regrid_controls.F90 \
    B16_trim.F90 \
    B17_checksum.F90 \
    B18_modele.F90
do

    echo "============================== Sampling $stage"
    xent -p $stage
done
