# Main program to plot all results from EntGVSD
# By: Nancy Kiang

#EntGVSD_ef.R
source("utils.R")

#***************************************************************************************************
# Step 2:  EDIT THIS SECTION FOR YOUR PATHS


# Choose data set version 
version = "v1.1"
entlclaidir = "lc_lai_ent/"
pathplot=paste(entlclaidir, "/plots/", sep="")

# Create output directory if it doesn't already exist
dir.create(pathplot)

# Step 3: Pick year
datatime=2004

#***************************************************************************************************
# Step 4:  Do plots. 
#   If inside R gui:  source("EntGVSD_ef.R")
#   If linux command line:  R CMD BATCH EntGVSD_ef.R 

# ent17
if (TRUE) {
    trimopt = "ent17"
    res="qxq"
    IM=1440
    JM=720
    filepre = "V1km_EntGVSD17G_BNUM"
    filepre = "V1km_EntGVSD17M_BNUM"
    filesuf = "_qxq"
    enttyp = 1:20
    map.entgvsd.steps(entlclaidir, res=res, enttyp,varname="lc", trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="laimax", trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    map.entgvsd.steps(entlclaidir, res=res,  enttyp, varname="hgt", trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    for (d in c( "2004_017", "2004_201")) {
        map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    }
}


# pure
if (TRUE) {
    trimopt = "pure"
    res="qxq"
    filepre = "V1km_EntGVSD16G_BNUM"
    filesuf = "_qxq"
    enttyp = 1:18
    map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="lc",  trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="laimax",  trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="hgt",  trimopt=trimopt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    for (d in c( "2004_Jan", "2004_Jul")) {
        map.entgvsd.steps(entlclaidir, res=res, enttyp=enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version=version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    }
}

# trimmming
if (TRUE) {
    trimopt = c("trimmed", "trimmed_scaled", "trimmed_scaled_nocrops", "trimmed_scaled_crops_ext")
    res="HXH"
    filepre = "VHXH_EntGVSD16G_BNUM"
    filesuf = "_plot"
    enttyp = 1:18
    for (opt in trimopt) {
        map.entgvsd.steps(entlclaidir, res=res, varname="lc", enttyp=enttyp, trimopt=opt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
        map.entgvsd.steps(entlclaidir, res=res, varname="laimax", enttyp=enttyp, trimopt=opt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
        map.entgvsd.steps(entlclaidir, res=res, varname="hgt",enttyp=enttyp, trimopt=opt, filepre, datatime, version, filesuf=filesuf,do.pdf = TRUE, pathplot=pathplot)
    }
}






