#!/bin/bash

exe=./merge.exe

#for year in {2003,2005,2006,2008,2009,2011,2012,2013}; do
year=2007

dir=/neponset/nbdata05/albedo/qsun/for_Nancy/2018/MCD12Q1/${year}
out_dir=/neponset/nbdata05/albedo/qsun/for_Nancy/2018/MCD12Q1/merged/${year}

if [ ! -r $out_dir ]; then
	mkdir -p $out_dir
fi

dsets=(0 4)
names=("IGBP" "PFT")

for ((i=0; i<2; i++)); do
	dset_i=${dsets[$i]}
	name=${names[$i]}

	list=""
	count=0

	for file in ${dir}/MCD12Q1*.hdf; do
		list="${list} $file"
		count=$(($count + 1))
	done

	if [ $count -lt 1 ]; then
		echo "No file found."
		exit
	fi

	out=${out_dir}/${name}_500m_${year}.bin

	echo "dset_i=$dset_i, name=$name, $count=$count, out=$out"
	$exe ${out} $dset_i $count $list

	if [ $? -ne 0 ]; then
		exit
	fi

	### reprojection ###
	cp out_500m_byte.hdr "${out}.hdr"
	out2=${out_dir}/${name}_500m_${year}_geo.bin
	gdalwarp -of ENVI -co INTERLEAVE=BIP -s_srs '+proj=sinu +a=6371007.181 +b=6371007.181 +units=m' -t_srs 'EPSG:4326' -r near -overwrite "${out}" "${out2}"

	gzip $out
	gzip $out2
done
#done






