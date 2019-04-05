#!/bin/bash
#
set -e








for stage in \
        A00_LAI3g_modis_entpftrevcrop.F90 \
        A01h_veg_height.f90 \
        A01_lc_laimax.f90 \
        A02_lc_lai_doy.f90 \
        A03_lc_lai_monthly.f90 \
        A04_trim_laimax_1kmx1km.F90 \
        A05_trim_laidoy_1kmx1km.F90 \
        A06_trim_lc_laimax_laimonth.F90 \
        A07_regrid.F90 \
        A08_trim.F90
do

    echo "============================== Running $stage"
    xent $stage
done

echo "=================== Finished run_all.sh"

