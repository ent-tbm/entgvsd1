#!/bin/sh

exe=./EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_subgrid_1km.exe
dir=/neponset/nbdata05/albedo/qsun/for_Nancy/2018/MCD12Q1/merged
dir_out=/neponset/nbdata05/albedo/qsun/for_Nancy/2018/MCD12Q1/merged/subgrid

#for year in {2002,2003,2005,2006,2008,2009,2011,2012,2013}; do
year=2007


	part=${dir}/${year}/PART_500m_${year}_geo.bin
	partz=${part}.gz
	if [ -r $partz ]; then
		echo "unzipping $part ..."
		gzip -d ${part}.gz
	fi
	if [ ! -r $part ]; then
		echo "$part not exists"
		exit
	fi

	d_out=${dir_out}/${year}
	if [ ! -r $d_out ]; then
		mkdir -p $d_out
	fi
	out=${d_out}/PART_SUB_1km_${year}_geo.hdf

	time $exe $part $out

	ret=$?
	gzip $part

	if [ $ret -ne 0 ]; then
		echo "$year failed"
		exit
	fi

#done
