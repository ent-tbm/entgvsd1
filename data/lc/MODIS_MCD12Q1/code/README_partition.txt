README:
Ent GVSD scripts and codes to partition MODIS IGBP and PFT land cover maps into maps of subgrid cover fractions of land cover types compatible with the Ent Terrestrial Biosphere Model.
 
Author: Qingsong Sun, Crystal Schaaf
 
DATA SOURCE:
MODIS MCD12Q1 IGBP and PFT land cover products at 500 m sinusoidal.
https://lpdaac.usgs.gov/products/mcd12q1v006/
 
Requirements:
Linux or unix system
gcc: 6.3.0
NetCDF: 4.3.0
HDF: 4.2.9
HDFEOS: 2.18
 
Contents:
1_merge.sh:  Run merge.c to merge MODIS tiled data into global data files.
 
2_partition.sh: Run EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_500m.c to partition the MODIS IGBP and PFT land cover into 29 cover classes compatible with Ent cover types.
 
3_subgrid.sh:   Run EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_subgrid_1km.c to organize the 500 m partitioning into 1 km subgrid fractions of the 29 cover classes.
 
EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_500m.c
 
EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_subgrid_1km.c
 
Makefile:  Compile the c programs.
 
merge.c:  Merge MODIS products in tiles to a global data file.
 
out_500m_byte.hdr:  ENVI metadata template.  Used by 1_merge.sh
 
HOW-TO:
	1	Install required software and libraries above.
	2	In your .profile or .bashrc, set environment variables for libraries:  df, netcdf, hdfeos.
	3	Get the MODIS MCD12Q1 IGBP and PFT files.
	4	Set up your output directory.
	5	In 1_merge.sh, edit the year(s) of data you wish to process.
	6	In 1_merge.sh, edit the dir for the path of the MODIS files.
	7	In 1_merge.sh, edit the out_dir for the path of your output directory.
	8	In 2_partition.sh, edit the dir to same path as your out_dir in 1_merge.sh.
***Actually, Qingsong, it might be nice to have 2_partition.sh have a different output directory from its input***
	9	In 2_partition.sh, edit the year(s) that you wish to process.
	10	mkdir your_subgrid_dir
	11	In 3_subgrid.sh, edit the dir to be the same as the output dir in 2_partition.sh.
	12	In 3_subgrid.sh, edit the dir_out to be your_subgrid_dir.
	13	make
	14	./1_merge.sh
	15	./2_partition.sh
	16	./3_subgrid.sh


