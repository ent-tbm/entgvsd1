#EntGVSD_efnk.R
source("utils_nk.R")

#***************************************************************************************************
# Step 1:  Add a subdirectory named PLOTS to the lc_lai_ent directory (so above the ent17, pure, etc).

#***************************************************************************************************
# Step 2:  EDIT THIS SECTION FOR YOUR PATHS

# Choose data set version 
version = "v1.1"

#Set paths:
#For Nancy on Mac:
if (FALSE) {
	path = "/Users/nkiang/NancyResearch/GISS/Models/Ent/Vegcover/Elizabeth_PLOTS/"
	#Choose date of output
	#dateout = "2019-07-19"
	#dateout = "2019-08-22"
	#dateout = "2019-08-26"
	#dateout = "2019-08-29"
	dateout = "2019-09-03"
	#entlclaidir = paste(dateout,"/lc_lai_ent.v7/", sep="")
	entlclaidir = paste(path, dateout,"/lc_lai_ent/", sep="")
        #If in R gui:
	#pathplot=paste(entlclaidir, "/PLOTS/", sep="")
        #If in git clone:
	pathplot="PLOTS/"
        dir.create(pathplot)
}
#For Nancy on gibbs:
if (FALSE) {
	entlclaidir = "/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent/"
	pathplot = "PLOTS/"	
	# Create output directory if it doesn't already exist
	dir.create(pathplot)
}
	
	
#For Elizabeth on gibbs:
if (TRUE) {
	entlclaidir = "lc_lai_ent/"
	pathplot=paste(entlclaidir, "/plots/", sep="")
	# Create output directory if it doesn't already exist
	dir.create(pathplot)
}

# Step 3: Pick year
datatime=2004

#***************************************************************************************************
# Step 4:  Select output to screen or to pdf.  Select whether to output checksums.
do.pdf = TRUE
do.checksum = TRUE
#***************************************************************************************************
# Step 5:  Do plots.
#
#   If inside R gui:  source("EntGVSD_ef.R")
#   If linux command line:  R CMD BATCH EntGVSD_ef.R 

# ent17
trimopt = "ent17"
res="qxq"
IM=1440
JM=720
filepre = "V1km_EntGVSD17G_BNUM"
filepre = "V1km_EntGVSD17M_BNUM"
filesuf = "_qxq"
enttyp = 1:20
map.entgvsd.steps(entlclaidir, res=res, enttyp,varname="lc", trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="laimax", trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
map.entgvsd.steps(entlclaidir, res=res,  enttyp, varname="hgt", trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
for (d in c( "2004_017", "2004_201")) {
	map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version, filesuf=filesuf,do.pdf=do.pdf, pathplot=pathplot, do.checksum = do.checksum)
}


# pure
trimopt = "pure"
res="qxq"
filepre = "V1km_EntGVSD16G_BNUM"
filesuf = "_qxq"
enttyp = 1:18
map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="lc",  trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="laimax",  trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = FALSE)
map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="hgt",  trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
for (d in c( "2004_Jan", "2004_Jul")) {
	map.entgvsd.steps(entlclaidir, res=res, enttyp=enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version=version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
}

# trimmming
trimopt = c("trimmed", "trimmed_scaled", "trimmed_scaled_nocrops", "trimmed_scaled_crops_ext")
res="HXH"
filepre = "VHXH_EntGVSD16G_BNUM"
filesuf = "_plot"
enttyp = 1:18
for (opt in trimopt) {
	map.entgvsd.steps(entlclaidir, res=res, varname="lc", enttyp=enttyp, trimopt=opt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
	map.entgvsd.steps(entlclaidir, res=res, varname="laimax", enttyp=enttyp, trimopt=opt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
	map.entgvsd.steps(entlclaidir, res=res, varname="hgt",enttyp=enttyp, trimopt=opt, filepre, datatime, version, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
}






