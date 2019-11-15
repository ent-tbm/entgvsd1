include mk/A00b_regrid_hr.mk
include mk/A00b_regrid_lr.mk
include mk/A00_LAI3g_modis_entpftrevcrop.mk
include mk/A01a_carrer_mean.mk
include mk/A01b_soil_albedo_hr.mk
include mk/A01b_soil_albedo.mk
include mk/A01f_albmodis_gridfill.mk
include mk/A01h_veg_heights.mk
include mk/A01_lc_laimax.mk
include mk/A02_lc_lai_doy.mk
include mk/A03_lc_lai_monthly.mk
include mk/A04_reclass_annual.mk
include mk/A05_reclass_doy.mk
include mk/A06_reclass_monthly.mk
include mk/A07_regrid_hr.mk
include mk/A07_regrid_lr.mk
include mk/A08_mc.mk
include mk/A08_nc.mk
include mk/A08_pu.mk
include mk/A08_tr.mk
include mk/A08_ts.mk


$(A00_LAI3g_modis_entpftrevcrop_OUTPUTS) : $(A00_LAI3g_modis_entpftrevcrop_INPUTS)
	xent -d A00_LAI3g_modis_entpftrevcrop.F90

$(A00b_regrid_hr_OUTPUTS) $(A00b_regrid_lr_OUTPUTS) : $(A00b_regrid_hr_INPUTS) $(A00b_regrid_hr_INPUTS)
    xent -d A00b_regrid.F90

$(A01a_carrer_mean_OUTPUTS) : $(A01a_carrer_mean_INPUTS)
	xent -d A01a_carrer_mean.F90

$(A01b_soil_albedo_hr_OUTPUTS) $(A01b_soil_albedo_hr_OUTPUTS) : $(A01b_soil_albedo__INPUTS) $(A01b_soil_albedo__INPUTS)
	xent -d A01b_soil_albedo.F90

$(A01f_albmodis_gridfill_OUTPUTS) : $(A01f_albmodis_gridfill_INPUTS)
	xent -d A01f_albmodis_gridfill.F90

$(A01h_veg_height_OUTPUTS) : $(A01h_veg_height_INPUTS)
	xent -d A01h_veg_height.F90

$(A01_lc_laimax_OUTPUTS) : $(A01_lc_laimax_INPUTS)
	xent -d A01_lc_laimax.F90

$(A02_lc_lai_doy_OUTPUTS) : $(A02_lc_lai_doy_INPUTS)
	xent -d A02_lc_lai_doy.F90

$(A03_lc_lai_monthly_OUTPUTS) : $(A03_lc_lai_monthly_INPUTS)
	xent -d A03_lc_lai_monthly.F90

$(A04_reclass_annual_OUTPUTS) : $(A04_reclass_annual_INPUTS)
	xent -d A04_reclass_annual.F90

$(A05_reclass_doy_OUTPUTS) : $(A05_reclass_doy_INPUTS)
	xent -d A05_reclass_doy.F90

$(A06_reclass_monthly_OUTPUTS) : $(A06_reclass_monthly_INPUTS)
	xent -d A06_reclass_monthly.F90

$(A07_regrid_hr_OUTPUTS) $(A07_regrid_lr_OUTPUTS) : $(A07_regrid_hr_INPUTS) $(A07_regrid_lr_INPUTS)
	xent -d A07_regrid.F90

$(A08_mc_OUTPUTS) $(A08_nc_OUTPUTS) $(A08_pu_OUTPUTS) $(A08_tr_OUTPUTS) $(A08_ts_OUTPUTS) : $(A08_mc_INPUTS) $(A08_nc_INPUTS) $(A08_pu_INPUTS) $(A08_tr_INPUTS) $(A08_ts_INPUTS) 
	xent -d A08_trim.F90
