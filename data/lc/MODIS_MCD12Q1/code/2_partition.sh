#!/bin/bash

exe=./EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_500m.exe
dir=/neponset/nbdata05/albedo/qsun/for_Nancy/2018/MCD12Q1/merged

#for year in {2003,2005,2006,2008,2009,2011,2012,2013}; do
year=2007

echo "PROCESSING $year ...";

for proj in {"","_geo"}; do

	IGBP=${dir}/${year}/IGBP_500m_${year}${proj}.bin
	PFT=${dir}/${year}/PFT_500m_${year}${proj}.bin
	out=${dir}/${year}/PART_500m_${year}${proj}.bin

	IGBPZ=${IGBP}.gz
	PFTZ=${PFT}.gz

	if [ -e $IGBPZ ]; then
		gzip -d $IGBPZ
	fi
	if [ -e $PFTZ ]; then
		gzip -d $PFTZ
	fi

	if [ ! -e $IGBP ] || [ ! -e $PFT ]; then
		echo "ERROR: IGBP=[$IGBP], PFT=[$PFT]"
		exit 1
	fi

	$exe $IGBP $PFT $out
	if [ $? -ne 0 ]; then
		echo "ERROR"
		exit 1
	fi

	hdr=${IGBP}.hdr
	hdr_out="${out%.bin}.hdr"
	if [ ! -e $hdr ]; then
		hdr=${IGBP%.bin}.hdr
		if [ ! -e $hdr ]; then
			echo "HDR NOT EXISTS: [$hdr]"
			exit 1
		fi
	fi

	cp $hdr ${hdr_out}

	echo "SUCCEED pft=[$PFT], igbp=[$IGBP], out=[$out], hdr_out=[$hdr_out]"

done

echo "ZIPPING $year ..."
gzip ${dir}/${year}/*.bin

#done
