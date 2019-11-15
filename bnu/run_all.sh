#!/bin/bash
#
set -e

# Stuck on B09...

B01_bnu_laimax.F90 \
B02_lc_laimax_modis_entpftrevcrop.F90 \
B03_regrid.F90 \
B04_veg_height.F90 \
B05_carrer_mean.F90 \
B06_albmodis_gridfill.F90 \
B07_soil_albedo.F90 \
B08_lc_laimax.F90 \
B09_lc_lai_doy.F90
##        A00a_bnu_laimax.F90 \
##        A00_LAI3g_modis_entpftrevcrop.F90 \
##        A00b_regrid.F90 \
##        A01h_veg_height.F90 \
##        A01a_carrer_mean.F90 \
##        A01f_albmodis_gridfill.F90 \
##        A01b_soil_albedo.F90 \
##        A01_lc_laimax.F90 \
##        A02_lc_lai_doy.F90 \
#        A03_lc_lai_monthly.F90 \
#        A04_reclass_annual.F90 \
#        A05_reclass_doy.F90 \
#        A06_reclass_monthly.F90 \


for stage in \
        A07_regrid.F90 \
        A07a_regrid_controls.F90 \
        A08_trim.F90 \
        A09_checksums.F90 \
        A10_modele.F90
do

    echo "============================== Running $stage"
    xent $stage
done

python3 A11_to_modele_format.py

# Plot everything
Rscript EntGVSD_ef.R

# Convert it all to PNG
python3 plots_to_png.py

echo "=================== Finished run_all.sh"

