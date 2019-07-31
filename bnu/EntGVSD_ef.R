#EntGVSD_ef.R
source("utils.R")
options(error=traceback)

print("AA1")

# Set paths and spatial resolution ----------------------------------------------------------------

#pathef = "/Users/nkiang/NancyResearch/GISS/Models/Ent/Vegcover/Elizabeth_PLOTS/"
pathef = "./"
res = "qxq"
IM=1440
JM=720

path = pathef

print("AA2")
# Choose data set version --------------------------------------------------------------------------
if (TRUE) { #Do Ent 17 PFTs
#Ent 17 PFTs
#dateout = "2019-07-19"
version = "v1.1"

entlclaidir = paste("lc_lai_ent/ent17/", sep="")
filepre = "V1km_EntGVSD17M_BNUM"
filepreht = filepre
datatime = "2004"
enttyp = 1:20  #Do all 20 layers.  paste(EntGVSD_COVER20, "          ")  #Elizabeth layers are 23 characters long.
if.trim = FALSE
trimopt = "raw"  #Please change this to "ent17"
pathplot=paste(pathef, dateout, "/PLOTS/", sep="")

} else {  # Do Ent 16 PFTs
	
---
#Ent 16 PFTs
#TBD
dateout = "2019-07-19"
version = "v1.1"

entlclaidir = paste(dateout,"/lc_lai_ent/", sep="")

}

print("AA3")
pathplot=paste(path, dateout, "/lc_lai_ent/PLOTS/", sep="")


#-------------------------------------------------------------------------------------------------
#1- Plot lc_max, lai_max,
##quartz(width=11,height=7) #Open plotting window
#quartz(width=11,height=6) #Open plotting window
zlim = c(0,1)
for (opt in trimopt) {
	if (if.trim) { #TRUE only for Ent 16 PFTs
		fname = paste(filepre, "_lc_max_", opt, "_", version, ".nc", sep="")
		pdf(file=paste(pathplot, fname, ".pdf", sep=""), width=11, height=7)
	} else {  #Ent 17 PFTs
		varname = "lc"
		fname = paste(filepre, "_",varname,"_",datatime,"_ann_raw_",version, "_", res,".nc", sep="")
	}
	filelc = paste(path, entlclaidir, fname, sep="")
	#par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
	par(mfrow=c(4,5), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
	titletop = paste(pathef, entlclaidir, sep="")
	map.EntGVSD.v1.1(filelc=filelc, file=filelc, res=res, zlim=zlim, varname=varname, layersnum=enttyp)
	entgvsd_pagetitles(titletop, fname)
	if (if.trim) {
		dev.off()
	}
}

for (opt in trimopt) {
	if (if.trim) { #TRUE only for Ent 16 PFTs
		fname = paste(filepre, "_lai_max_", opt, "_", version, ".nc", sep="")
		pdf(file=paste(pathplot, fname, ".pdf", sep=""), width=11, height=7)
	} else {
		varname = "laimax"
		fname = paste(filepre, "_",varname,"_",datatime,"_raw_",version, "_", res,".nc", sep="")	}
	zlim = c(0,7)
	filein = paste(path, entlclaidir, fname, sep="")
	print(filein)
	#par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
	par(mfrow=c(4,5), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
	map.EntGVSD.v1.1(filelc=filelc, file=filein, varname=varname, res=res, zlim=zlim, layersnum=enttyp, colors=drywet(40))
	entgvsd_pagetitles(paste(pathef, entlclaidir, sep=""), fname)
	if (if.trim) {
		dev.off()
	}
}

#2- Plot "trimmed_scaled_crops_ext1"
par(mfrow=c(2,2), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
if (if.trim) { #TRUE only for Ent 16 PFTs
opt = "trimmed_scaled_crops_ext1"
for (varname in c("lc", "laimax", "hgt")) {
	fname = paste(filepre, "_", varname, "_2004_", opt, "_", version, "_", res,".nc", sep="")
	#pdf(file=paste(pathplot, fname, ".pdf", sep=""), width=11, height=7)
	zlim = c(0,7)
	filein = paste(path, entlclaidir, fname, sep="")
	print(filein)
	#par(mfrow=c(2,2), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
	ncid = open.nc(con=filein, write=FALSE)
	x = var.get.nc(ncid, paste(varpre,"crops_herb", sep=""))
	layers = var.get.nc(ncid, "layers")
	x[x==0]=NA
	plot.grid.continuous(mapz=x[,,match("crops_herb", trim(layers))], zlim=zlim, res=res, colors=drywet(40))
		plot(coastsCoarse, add=TRUE, lwd=0.25)
		title(paste("  ", varname, "crops_herb ext"))
	cex = 1
	text(-175, 0, varname, adj=0, cex=cex)
	text(-175, -20, paste("mean", round(na.mean(x),2)), adj=0, cex=cex)
	text(-175, -35, paste("max", round(na.max(x),2)), adj=0, cex=cex)
	#map.EntGVSD(filelc=filelc, file=filein, res=res, zlim=zlim, varlist="crops_herb", 	colors=drywet(40))
	#dev.off()
}}
entgvsd_pagetitles(paste(pathef, entlclaidir, sep=""), filepre)
	

#3- Plot height 
for (opt in trimopt) {
	varname = "hgt"
	if (if.trim) { #Ent 16 PFTs
	fname = paste(filepre, "_", varname, "_2004_", opt, "_", version, "_", res,".nc", sep="")
		pdf(file=paste(pathplot, fname, ".pdf", sep=""), width=11, height=7)
		varpre=""
	} else {
		fname = paste(filepre, "_", varname, "_2004_", opt, "_", version, "_", res,".nc", sep="")
	}
	filein = paste(path, entlclaidir, fname, sep="")
	print(filein)
	#par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
	par(mfrow=c(4,5), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
	zlim = c(0,40)
	map.EntGVSD.v1.1(filelc=filelc, file=filein, res=res, zlim=zlim, varname=varname, 	layersnum=enttyp)
	entgvsd_pagetitles(path, fname)
	if (if.trim) {
		dev.off()
	}
}

#4- Plot lai on day of year
for (opt in trimopt) {
	for (doy in c("017", "201")) {
	if (if.trim) { #Ent 16 PFTs
		fname = paste(filepre, "_lai_2004_",doy,"_", opt, "_", version, "_", res,".nc", sep="")
		pdf(file=paste(pathplot, fname, ".pdf", sep=""), width=11, height=6)
	} else { #Ent 17 PFTs
		fname = paste(filepre, "_lai_2004_", doy,"_",opt,"_",version,"_", res,".nc", sep="")
		pdf(file=paste(pathplot, fname,".pdf", sep=""), width=11, height=6)
	}
	par(mfrow=c(4,5), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)

		filein = paste(path, entlclaidir, "ent17/", fname, sep="")	
		print(filein)
		zlim=c(0,7)
		map.EntGVSD.v1.1(filelc=filelc, file=filein, varname="lai", layersnum=enttyp, res=res, zlim=zlim)
		entgvsd_pagetitles(paste(path, entlclaidir, fname, sep=""))
	close.nc(ncid)
	dev.off()
	}
}



