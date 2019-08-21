#!/bin/bash
#
set -e

#        A00_LAI3g_modis_entpftrevcrop.F90 \
#        A01h_veg_height.F90 \
#        A01_lc_laimax.F90 \
#        A02_lc_lai_doy.F90 \
#        A03_lc_lai_monthly.F90 \

for stage in \
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

