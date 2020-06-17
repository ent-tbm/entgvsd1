#!/bin/bash

version=v1.1
year=2004
timestep=annual
#Select whether to plot filled of no-fill soil albedo
#gap-filled
# ifill=_fill
# ifilldir=_fill
#not gap-filled
 ifill=''
 ifilldir=_nofill

path=${PWD}
DIR="$(dirname $path)"/outputs/soilalbedo$ifilldir
echo $DIR

RES[1]=2HX2
RES[2]=HXH
RES[3]=6km
RES[4]=1km

NRES=3
k=0
while [ $k -lt $NRES ]; do
let k=k+1
Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_bs_brightratio$ifill.nc bs_brightratio
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_2HX2_bs_brightratio$ifill.nc bs_brightratio
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_bs_brightratio$ifill.nc bs_brightratio
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_6km_bs_brightratio$ifill.nc bs_brightratio
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_1km_bs_brightratio$ifill.nc bs_brightratio

Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_EntGVSD_${version}_CarrerGISS_VIS_${timestep}_${year}$ifill.nc soilalbedo_VIS
Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_EntGVSD_${version}_CarrerGISS_NIR1_${timestep}_${year}$ifill.nc soilalbedo_NIR1
Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_EntGVSD_${version}_CarrerGISS_NIR2_${timestep}_${year}$ifill.nc soilalbedo_NIR2
Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_EntGVSD_${version}_CarrerGISS_NIR3_${timestep}_${year}$ifill.nc soilalbedo_NIR3
Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_EntGVSD_${version}_CarrerGISS_NIR4_${timestep}_${year}$ifill.nc soilalbedo_NIR4
Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_EntGVSD_${version}_CarrerGISS_NIR5_${timestep}_${year}$ifill.nc soilalbedo_NIR5
Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_${RES[k]}_EntGVSD_${version}_CarrerGISS_SW_${timestep}_${year}$ifill.nc soilalbedo_SW

#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_EntGVSD_${version}_CarrerGISS_VIS_${timestep}_${year}$ifill.nc soilalbedo_VIS
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_EntGVSD_${version}_CarrerGISS_NIR1_${timestep}_${year}$ifill.nc soilalbedo_NIR1
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_EntGVSD_${version}_CarrerGISS_NIR2_${timestep}_${year}$ifill.nc soilalbedo_NIR2
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_EntGVSD_${version}_CarrerGISS_NIR3_${timestep}_${year}$ifill.nc soilalbedo_NIR3
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_EntGVSD_${version}_CarrerGISS_NIR4_${timestep}_${year}$ifill.nc soilalbedo_NIR4
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_EntGVSD_${version}_CarrerGISS_NIR5_${timestep}_${year}$ifill.nc soilalbedo_NIR5
#Rscript B22_plots_custom.R 16 soilalbedo $DIR soilalbedo_HXH_EntGVSD_${version}_CarrerGISS_SW_${timestep}_${year}$ifill.nc soilalbedo_SW

#Rscript B22_plots_custom.R 16 lc $DIR/pure V1km_EntGVSD_${version}_16G_BMSa_lc_${year}_ann_pure$ifill.nc

done

