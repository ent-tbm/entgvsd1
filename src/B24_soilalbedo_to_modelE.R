#B24_soilalbedo_to_modelE.R
# Reformats netcdf individual month LAI files into a single file with lai(im, jm, lc, month)
# To run: R CMD BATCH B24_soilalbedo_to_modelE.R

source("utils.R")
library("RNetCDF")

#Update 
FILL = '_fill'; FILLdescription="  FILLED version: Undefined cells interpolated."
#FILL = ''; FILLdescription="  No interpolation."
YEAR='2004'
version = 'v1.1'

reschar = c("2HX2", "HXH")
#res="2HX2"
#res="HXH"

for (res in reschar) {

path = '../outputs/soilalbedo/'
file = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_multiband_annual_',YEAR,FILL,'.nc', sep='')
fileout = paste(path, file, sep='')
print(fileout)

description = paste('Ent Global Vegetation Structure Dataset (Ent GVSD) Soil albedo derived from D. Carrer et al. (2014) MODIS soil albedos.  Created with Ent GVSD ',version,'.',FILLdescription,  sep="")

create.soilalbedo.ModelE.template.nc(res=res, description, fileout, contact="Nancy.Y.Kiang@nasa.gov")

nco = open.nc(con=fileout, write=T)

#lon, lat
filein = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_SW_annual_',YEAR,FILL,'.nc', sep='')
nci = open.nc(con=paste(path,filein,sep=""))

var.put.nc(nco, 'lon', var.get.nc(nci,'lon'))
var.put.nc(nco, 'lat', var.get.nc(nci,'lat'))

# SW
filein = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_SW_annual_',YEAR,FILL,'.nc', sep='')
print(filein)
nci = open.nc(con=paste(path,filein,sep=""))
#var.put.nc(nco, 'albedo_soil_SW', var.get.nc(nci, 'albSW'))
var.put.nc(nco, 'albedo_soil_SW', var.get.nc(nci, 'soilalbedo_SW'))
close.nc(nci)

# GISS Bands
GISSBAND = var.get.nc(nco, 'band_names')
dim3D = dim(var.get.nc(nco, 'albedo_soil_giss'))
print(dim3D)
alb3D = array(NA, dim=dim3D)
for (i in 1:length(GISSBAND)) {
   filein = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_',trim(GISSBAND[i]),'_annual_',YEAR,FILL,'.nc', sep='')
   print(filein)
   nci = open.nc(con=paste(path,filein,sep=""))
#   varin = paste('albgiss_', trim(GISSBAND[i]),sep='')
   varin = paste('soilalbedo_', trim(GISSBAND[i]),sep='')
   print(varin)
   alb3D[,,i] = var.get.nc(nci, varin)  
   
   close.nc(nci)
}
 
   var.put.nc(nco, 'albedo_soil_giss', alb3D)
   close.nc(nco)

} #for res


