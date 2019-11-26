#! /bin/bash

# 1) Interpolate monthly CRU temperature from 0.5º x 0.5º to 1km x 1km

cd ./CRU_GPCC

ulimit -s unlimited

module purge
module load other/comp/gcc-4.9.2-sp3
module load other/ncl-6.3.0

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include CRU_Temp_Interp_to_1kmx1km_monthly.f

gfortran -o myExe arrayutil.o CRU_Temp_Interp_to_1kmx1km_monthly.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

./myExe


# 2) Interpolate monthly GPCC precipitation from 0.5º x 0.5º to 1km x 1km

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include GPCC_Precip_Interp_to_1kmx1km_monthly.f

gfortran -o myExe arrayutil.o GPCC_Precip_Interp_to_1kmx1km_monthly.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

./myExe

cd ..


# 3) Calculate average of climate metrics for C4 classification

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include CRU_GPCC_C4.f

gfortran -o myExe arrayutil.o CRU_GPCC_C4.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

./myExe


# 4) Climate classification based on gridded CRU temperature and GPCC precipitation

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include arrayutil.f

gfortran -c -cpp -mcmodel=medium -fconvert=big-endian -O2 -fno-range-check -I/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/include CRU_GPCC_c4norm.f

gfortran -o myExe arrayutil.o CRU_GPCC_c4norm.o -L/usr/local/other/SLES11.3/netcdf/3.6.3/gcc-4.9.2-sp3/lib -lnetcdf

./myExe

