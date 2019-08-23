#!/bin/bash
#
set -e


for stage in \
        A00_LAI3g_modis_entpftrevcrop.F90 \
        A00b_regrid.F90 \
        A01h_veg_height.F90 \
        A01a_carrer_mean.F90 \
        A01b_soil_albedo.F90 \
        A01_lc_laimax.F90 \
        A02_lc_lai_doy.F90 \
        A03_lc_lai_monthly.F90 \
        A04_reclass_annual.F90 \
        A05_reclass_doy.F90 \
        A06_reclass_monthly.F90 \
        A07_regrid.F90 \
        A08_trim.F90
do

    echo "============================== Running $stage"
    xent $stage
done

echo "=================== Finished run_all.sh"

