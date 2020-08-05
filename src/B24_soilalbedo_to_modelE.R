#B24_soilalbedo_to_modelE.R
# Reformats netcdf individual band soil albedo files into a single file with albedo_soil_giss(im, jm, band) and albedo_soil_SW(im, jm)
# To run: Rscript B24_soilalbedo_to_modelE.R <fill | nofill>

source("utils.R")
library("RNetCDF")

args = commandArgs(trailingOnly=TRUE)
print(args)
numargs = length(args)
if (numargs != 1) {
  print ('Usage:  Rscript B24_soilalbedo_to_modelE.R <fill | nofill>')
  print ('Input fill or nofill to convert soil albedo files that are or are not gap-filled.')
  quit()
} 

if (args[1]=='fill') {
   #gap-filled
   IFILL = '_fill'; FILLdescription='  Gap-filled.'
   IFILLdir = '_fill'
} else {
   #no gap-filled
   IFILL = ''; FILLdescription='  Not gap-filled.'
   IFILLdir = '_nofill'
}

YEAR='2004'
version = 'v1.1'

reschar = c("2HX2", "HXH")
#res="2HX2"
#res="HXH"

for (res in reschar) {

path = paste('../outputs/soilalbedo',IFILLdir,'/', sep="")
file = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_multiband_annual_',YEAR,IFILL,'.nc', sep='')
fileout = paste(path, file, sep='')
print(fileout)
file2d = paste('soilalbedo_',res,'_EntGVSD_',version,'_Carrer_GISSbands_annual_',YEAR,IFILL,'.nc', sep='')
fileout2d = paste(path, file2d, sep='')
print(fileout2d)

description = paste('Ent Global Vegetation Structure Dataset (Ent GVSD) Soil albedo derived from D. Carrer et al. (2014) MODIS soil albedos.  Created with Ent GVSD ',version,'.',FILLdescription,  sep="")

create.soilalbedo.ModelE.template.nc(res=res, description, fileout, contact="Nancy.Y.Kiang@nasa.gov")
create.soilalbedo.ModelE.template.2D.nc(res, description, fileout2d, contact="Nancy.Y.Kiang@nasa.gov")

nco = open.nc(con=fileout, write=T)
nco2d = open.nc(con=fileout2d, write=T)

#lon, lat
filein = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_SW_annual_',YEAR,IFILL,'.nc', sep='')
nci = open.nc(con=paste(path,filein,sep=""))

var.put.nc(nco, 'lon', var.get.nc(nci,'lon'))
var.put.nc(nco, 'lat', var.get.nc(nci,'lat'))
var.put.nc(nco2d, 'lon', var.get.nc(nci,'lon'))
var.put.nc(nco2d, 'lat', var.get.nc(nci,'lat'))

# SW
filein = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_SW_annual_',YEAR,IFILL,'.nc', sep='')
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
   filein = paste('soilalbedo_',res,'_EntGVSD_',version,'_CarrerGISS_',trimws(GISSBAND[i]),'_annual_',YEAR,IFILL,'.nc', sep='')
   print(filein)
   nci = open.nc(con=paste(path,filein,sep=""))
#   varin = paste('albgiss_', trim(GISSBAND[i]),sep='')
   varin = paste('soilalbedo_', trimws(GISSBAND[i]),sep='')
   print(varin)
   alb3D[,,i] = var.get.nc(nci, varin)  
   var.put.nc(nco2d, paste('soilalb_',trimws(GISSBAND[i]), sep=""), alb3D[,,i]) 
   close.nc(nci)
}
 
   var.put.nc(nco, 'albedo_soil_giss', alb3D)
   close.nc(nco)
   close.nc(nco2d)

} #for res


