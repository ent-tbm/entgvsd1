#EntGVSD_efnk.R
source("utils.R")

#******************************************************************************************
# Step 1:  EDIT THIS SECTION FOR YOUR PATHS AND TO SET UP A plots DIRECTORY.

# Choose data set version 
version = "v1.1"

#Set paths:

#For Nancy on gibbs:
if (FALSE) {
    outputsdir = "/home2/rpfische/git/entgvsd1/outputs/"
    pathplot = "../outputs/plots/" 
    # Create output directory if it doesn't already exist
    dir.create("../outputs/")
    dir.create(pathplot)
}
    
    
#For running within the src/ directory as part of the generation process:
if (TRUE) {
    outputsdir = "../outputs/"
    pathplot=paste(outputsdir, "plots/", sep="")
    # Create output directory if it doesn't already exist
    dir.create(pathplot)
}

# Step 2: Pick year
datatime=2004

#******************************************************************************************
# Step 3:  Select output to screen or to pdf.  Select whether to output checksums.
do.pdf = TRUE
do.checksum = TRUE
#******************************************************************************************
# Step 4:  Do plots.
#
#   If inside R gui:  source("EntGVSD_ef.R")
#   If linux command line:  R CMD BATCH EntGVSD_ef.R 

# ent17
if (TRUE) {
    trimopt = "ent17"
    res="qxq"
    IM=1440
    JM=720
    #filepre = "V1km_EntGVSD17G_BNUM"
    filepre = "V1km_EntGVSD"
    icov = "17M"
    idat = "BMSa"
    
    filesuf = "_forplot"
    enttyp = 1:20

    map.entgvsd.steps(outputsdir, res=res, enttyp,varname="lc", trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(outputsdir, res=res, enttyp, varname="laimax", trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(outputsdir, res=res,  enttyp, varname="hgt", trimopt=trimopt, "", filepre, datatime, version,  icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    for (d in c( "2004_017", "2004_201")) {
        map.entgvsd.steps(outputsdir, res=res, enttyp, varname="lai", trimopt=trimopt, "", filepre, datatime=d, version, icov, idat, filesuf=filesuf,do.pdf=do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }

    for (m in 1:12) {
        d = paste(datatime, "_", MON[m], sep="")
        outpre = paste(MONnum[m],"_", sep="")
        map.entgvsd.steps(outputsdir, res=res, enttyp, varname="lai", trimopt=trimopt, outpre, filepre, datatime=d, version, icov, idat, filesuf=filesuf,do.pdf=do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }

	#Layered check (only for ent17):
    map.entgvsd.steps(outputsdir, res=res, enttyp, varname="laimax_err", trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = FALSE) #This is by PFT

    #Single layer checks
    # Monthly
    for (m in 1:12) {
        d = paste(datatime, "_", MON[m], sep="")
        outpre = paste(MONnum[m],"_", sep="")
    	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    	
    }
	
    # Totals
    # Dominant lc
    fname = paste(trimopt,"/", filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
    domlc = Ent_calc_domlc(file=paste(outputsdir, fname, sep=""), enttyp)
    fnameout = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc_domlc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
    pdf(file=paste(pathplot, fnameout, ".pdf", sep=""), width=8, height=5)
    Ent_domlc_plot(lctype=domlc, numpft=17, res=res, legend.cex=0.6, Entcolors=Entcolors17[match(na.min(domlc), Entcolors17[,"num"]):21,], if.new=FALSE)
    mtext(fnameout, cex=0.8)
    dev.off()

    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompft", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    #Redundant plot, just test if this function works, too.
    #map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #Colors do not work in R as categorical
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompftlc", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 

	#npftgrid
	fname = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	#file = paste(outputsdir, fname, sep="")
	fileout = paste(pathplot, filepre, "_",version, "_",icov, "_",idat, "_","npftgrid", "_",datatime,"_ann",  "_", trimopt, filesuf, ".pdf", sep="")
	pdf(file=fileout, width=8, height=5)
	npftgrid.list = Ent_calc_npftgrid(pathin=paste(outputsdir,trimopt,"/", sep=""), fname, pathout=pathplot, npft=16)
        title(npftgrid.list$file.nc)
	dev.off()
	
     map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.
	
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.

   #map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclaimax_err", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
    
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.
 
 	#map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lchgt_err", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
  	
  	#Checksum diff maps
  	res="HXH"
  	checksumdir = paste(outputsdir, "checksum/", sep="")
  	checksuff = paste("_diff", filesuf, sep="")
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 

	for (m in 1:12) {
           d = paste(datatime, "_", MON[m], sep="")
           outpre = paste(MONnum[m],"_", sep="")
	   map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
    for (d in c( "2004_017", "2004_201")) {
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, "", filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 

	
} #ent17

# pure
if (TRUE) {
    trimopt = "pure"
    res="qxq"
    filepre = "V1km_EntGVSD"
    icov = "16G"
    idat = "BMSa"
    filesuf = "_forplot"
    enttyp = 1:18
    
    map.entgvsd.steps(outputsdir, res=res,  enttyp=enttyp, varname="lc",  trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(outputsdir, res=res,  enttyp=enttyp, varname="laimax",  trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(outputsdir, res=res,  enttyp=enttyp, varname="hgt",  trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    #for (d in c( "2004_Jan", "2004_Jul")) {

    for (m in 1:12) {
        d = paste(datatime, "_", MON[m], sep="")
        outpre = paste(MONnum[m],"_", sep="")
        map.entgvsd.steps(outputsdir, res=res, enttyp=enttyp, varname="lai", trimopt=trimopt, outpre, filepre, datatime=d, version=version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }
 

    #Single layer checks
    # Monthly

    for (m in 1:12) {
        d = paste(datatime, "_", MON[m], sep="")
        outpre = paste(MONnum[m],"_", sep="")
    	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    
	}
		
	# Totals
	# Dominant lc
	fname = paste(trimopt,"/", filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	domlc = Ent_calc_domlc(file=paste(outputsdir, fname, sep=""), enttyp)
	fnameout = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc_domlc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	pdf(file=paste(pathplot, fnameout, ".pdf", sep=""), width=8, height=5)
	Ent_domlc_plot(lctype=domlc, numpft=16, res=res, legend.cex=0.6, Entcolors=Entcolors16[match(na.min(domlc), Entcolors16[,"num"]):dim(Entcolors16)[1],], if.new=FALSE)
	mtext(fnameout, cex=0.8)
	dev.off()
   	
	if (FALSE) {
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompft", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    	
        #map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompftlc", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	} #if (FALSE)
 	
	#npftgrid
	fname = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	#file = paste(outputsdir, fname, sep="")
	fileout = paste(pathplot, filepre, "_",version, "_",icov, "_",idat, "_","npftgrid", "_",datatime,"_ann",  "_", trimopt, filesuf, ".pdf", sep="")
	pdf(file=fileout, width=8, height=5)
	npftgrid.list = Ent_calc_npftgrid(pathin=paste(outputsdir,trimopt,"/", sep=""), fname, pathout=pathplot, npft=16)
        title(npftgrid.list$file.nc)
	dev.off()
	
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.

	
  	#Checksum diff maps
  	res = "HXH"
  	checksumdir = paste(outputsdir, "checksum/",  sep="")
  	checksuff = paste("_diff", filesuf, sep="")
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 

        for (m in 1:12) {
           d = paste(datatime, "_", MON[m], sep="")
           outpre = paste(MONnum[m],"_", sep="")
     	   map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
    for (d in c( "2004_017", "2004_201")) {
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, "", filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 

} #pure

 
# purelr
if (TRUE) {
    trimopt = "purelr"
    res="HXH"
    filepre = "VHXH_EntGVSD"
    icov = "16G"
    idat = "BMSa"
    filesuf = ""
    enttyp = 1:18     
  	do.checksum=FALSE

    map.entgvsd.steps(outputsdir, res=res,  enttyp=enttyp, varname="lc",  trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(outputsdir, res=res,  enttyp=enttyp, varname="laimax",  trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(outputsdir, res=res,  enttyp=enttyp, varname="hgt",  trimopt=trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    #for (d in c( "2004_Jan", "2004_Jul")) {
    for (m in 1:12) {
        d = paste(datatime, "_", MON[m], sep="")
        outpre = paste(MONnum[m],"_", sep="")
        map.entgvsd.steps(outputsdir, res=res, enttyp=enttyp, varname="lai", trimopt=trimopt, outpre, filepre, datatime=d, version=version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }
 
    #Single layer checks
    if (FALSE) { #----
    # Monthly
    for (m in 1:12) {
        d = paste(datatime, "_", MON[m], sep="")
        outpre = paste(MONnum[m],"_", sep="")
    	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    
	}
    } #if FALSE ----

   	
	if (FALSE) { #-------
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompft", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    	
        #map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompftlc", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	} #if (FALSE) ---------
 	
    # Totals
    # Dominant lc
    fname = paste(trimopt,"/", filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
    domlc = Ent_calc_domlc(file=paste(outputsdir, fname, sep=""), enttyp)
    fnameout = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc_domlc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
    pdf(file=paste(pathplot, fnameout, ".pdf", sep=""), width=8, height=5)
    Ent_domlc_plot(lctype=domlc, numpft=17, res=res, legend.cex=0.6, Entcolors=Entcolors17[match(na.min(domlc), Entcolors17[,"num"]):21,], if.new=FALSE)
    mtext(fnameout, cex=0.8)
    dev.off()

 	if (do.checksum) {
	#npftgrid
	fname = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	#file = paste(outputsdir, fname, sep="")
	fileout = paste(pathplot, filepre, "_",version, "_",icov, "_",idat, "_","npftgrid", "_",datatime,"_ann",  "_", trimopt, filesuf, ".pdf", sep="")
	pdf(file=fileout, width=8, height=5)
	npftgrid.list = Ent_calc_npftgrid(pathin=paste(outputsdir,trimopt,"/", sep=""), fname, pathout=pathplot, npft=16)
        title(npftgrid.list$file.nc)
	dev.off()
	
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.
	} #--------

	# Totals
	# Dominant lc
	fname = paste(trimopt,"/", filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	domlc = Ent_calc_domlc(file=paste(outputsdir, fname, sep=""), enttyp)
	fnameout = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc_domlc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	pdf(file=paste(pathplot, fnameout, ".pdf", sep=""), width=8, height=5)
	Ent_domlc_plot(lctype=domlc, numpft=16, res=res, legend.cex=0.6, Entcolors=Entcolors16[match(na.min(domlc), Entcolors16[,"num"]):dim(Entcolors16)[1],], if.new=FALSE)
	mtext(fnameout, cex=0.8)
	dev.off()
	
	#lc checksum in R
	fname= paste(filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	#file = paste(outputsdir, fname, sep="")
	lcchecksum.list = Ent_calc_lc_checksum(pathin=paste(outputsdir, trimopt,"/", sep=""), fname, pathout=pathplot, enttyp=enttyp)
	fileout = paste(pathplot, filepre, "_",version, "_",icov, "_",idat, "_","lc_checksum", "_",datatime, trimopt, filesuf, ".pdf", sep="")
	pdf(file=fileout, width=8, height=5)
    #map.GCM(file=paste(file,"_checksum.nc", sep=""), res=res, varname="lc_checksum")  
    map.GCM(file=lcchecksum.list$file.nc, res=res, varname="lc_checksum")  
    title(lcchecksum.list$file.nc)
	dev.off()
	
	#bs_brightratio
	fname = 'bs_brightratio.nc'
	file = paste(outputsdir, trimopt,"/",fname, sep="")
	fileout = paste(pathplot, filepre, "_",version, "_",icov, "_",idat, "_","bs_brightratio", "_",datatime, trimopt, filesuf, ".pdf", sep="")
	pdf(file=fileout, width=8, height=5)
    map.GCM(file=file, res=res, colors=giss.palette.nowhite(40), varname="bs_brightratio")  
    title(file)
	dev.off()
	
	#npftgrid
	fname = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	#file = paste(outputsdir, fname, sep="")
	fileout = paste(pathplot, filepre, "_",version, "_",icov, "_",idat, "_","npftgrid", "_",datatime,"_ann",  "_", trimopt, filesuf, ".pdf", sep="")
	pdf(file=fileout, width=8, height=5)
	npftgrid.list = Ent_calc_npftgrid(pathin=paste(outputsdir,trimopt,"/", sep=""), fname, pathout=pathplot, npft=16)
        title(npftgrid.list$file.nc)
	dev.off()
	

  	#Checksum diff maps
  	if (FALSE) { #---------------------
  	res = "HXH"
  	checksumdir = paste(outputsdir, "checksum/",  sep="")
  	checksuff = paste("_diff", filesuf, sep="")
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
        for (m in 1:12) {
           d = paste(datatime, "_", MON[m], sep="")
           outpre = paste(MONnum[m],"_", sep="")
	   map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
        for (d in c( "2004_017", "2004_201")) {
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, "", filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
  	 
  	 } # if FALSE -----------------

} #purelr


do.checksum = TRUE
# trimmming
if (TRUE) {
   trimopts = c("trimmed", "trimmed_scaled", "trimmed_scaled_nocrops", "trimmed_scaled_crops_ext")
    res="HXH"
    filepre = "VHXH_EntGVSD"
    icov = "16G"
    idat = "BMSa"
    filesuf = "_forplot"
    enttyp = 1:18
    
    for (trimopt in trimopts) {
    		
  	for (varname in c("lc", "laimax", "hgt")) {
	        map.entgvsd.steps(outputsdir, res=res, varname, enttyp=enttyp, trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
	        #map.entgvsd.steps(outputsdir, res=res, varname="lc", enttyp=enttyp, trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    	    #map.entgvsd.steps(outputsdir, res=res, varname="laimax", enttyp=enttyp, trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
        	#map.entgvsd.steps(outputsdir, res=res, varname="hgt",enttyp=enttyp, trimopt, "", filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
	}
    	
	for (m in 1:12) {
	  d = paste(datatime, "_", MON[m], sep="")
	  outpre = paste(MONnum[m],"_", sep="")
	  map.entgvsd.steps(outputsdir, res=res, enttyp=enttyp, varname="lai", trimopt, outpre, filepre, datatime=d, version=version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    	}
    
   		#Single layer checks
   	 	# Monthly
		for (m in 1:12) {
		   d = paste(datatime, "_", MON[m], sep="")
		   outpre = paste(MONnum[m],"_", sep="")
    		   map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    
		}

	    	
		# Totals
		fname = paste(trimopt,"/", filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
		domlc = Ent_calc_domlc(file=paste(outputsdir, fname, sep=""), enttyp)
		fnameout = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc_domlc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
		pdf(file=paste(pathplot, fnameout, ".pdf", sep=""), width=8, height=5)
		Ent_domlc_plot(lctype=domlc, numpft=16, res=res, legend.cex=0.6, Entcolors=Entcolors16[match(na.min(domlc), Entcolors16[,"num"]):dim(Entcolors16)[1],], if.new=FALSE)
		mtext(fnameout, cex=0.8)
		dev.off()

		#npftgrid
		fname = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
		fileout = paste(pathplot, filepre, "_",version, "_",icov, "_",idat, "_","npftgrid", "_",datatime,"_ann",  "_", trimopt, filesuf, ".pdf", sep="")
		pdf(file=fileout, width=8, height=5)
		npftgrid.list = Ent_calc_npftgrid(pathin=paste(outputsdir,trimopt,"/", sep=""), fname, pathout=pathplot, npft=16)
	    title(npftgrid.list$file.nc)
		dev.off()
	
		if (FALSE) {
   	 	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompft", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    
    	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 		map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_dompftlc", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 		} #if (FALSE)
 		
    	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lc_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    	map.entgvsd.check.misc(outputsdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.
 	
 	#Checksum diff maps
  	checksumdir = paste(outputsdir, "checksum/", sep="")
  	checksuff = "_diff"
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	for (m in 1:12) {
	   d = paste(datatime, "_", MON[m], sep="")
	   outpre = paste(MONnum[m],"_", sep="")
  	   map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, outpre, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
 	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, "", filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 		
   } #trimopts
} #trimming

#modelE
if (TRUE) {
    trimopt = c("modelE")
    res="2HX2"
    filepre = "V2HX2_EntGVSD"
    icov = "16G"
    idat = "BMSa"
    enttyp = 1:18
    covlist = EntGVSD_COV13

    filesuf = paste("_ann_",trimopt,sep="")

    if (do.pdf) {

        varname = "lc"
        fname = paste(filepre, "_", version, "_", icov, "_", idat, "_", varname,"_", datatime,filesuf, ".nc3", sep="")
        file = paste(outputsdir,trimopt,"/",fname, sep="")
        print(file)
        filepdf = paste(pathplot, fname, ".pdf", sep="")
        pdf(file=filepdf, width=11, height=7)
        par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
        map.EntGVSD(filelc=NULL, file, res=res, varpre="", varlist=EntGVSD_COV13, colors=giss.palette.nowhite(40), type="any", zlim=c(0,1), if.zeroNA=TRUE, titletype=1)
        mtext(outer=TRUE, file)
        dev.off()

        varname = "laimax"
        fname = paste(filepre, "_", version, "_", icov, "_", idat, "_", varname,"_", datatime,filesuf, ".nc3", sep="")
        file = paste(outputsdir,trimopt,"/",fname, sep="")
        print(file)
        filepdf = paste(pathplot, fname, ".pdf", sep="")
        pdf(file=filepdf, width=11, height=7)
        par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
        map.EntGVSD(filelc=NULL, file, res=res, varpre="", varlist=EntGVSD_COV13, colors=drywet(40), type="any", zlim=c(0,7), if.zeroNA=TRUE, titletype=1)
        mtext(outer=TRUE, file)
        dev.off()

        varname = "hgt"
        fname = paste(filepre, "_", version, "_", icov, "_", idat, "_", varname,"_", datatime,filesuf, ".nc3", sep="")
        file = paste(outputsdir,trimopt,"/",fname, sep="")
        print(file)
        filepdf = paste(pathplot, fname, ".pdf", sep="")
        pdf(file=filepdf, width=11, height=7)
        par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
        map.EntGVSD(filelc=NULL, file, res=res, varpre="hgt_", varlist=EntGVSD_COV13, colors=giss.palette.nowhite(40), type="any", zlim=c(0,40), if.zeroNA=TRUE, titletype=1) 
        mtext(outer=TRUE, file)
        dev.off()

        varname = "lai"
        filesuf=paste("_", trimopt, sep="")
        fname = paste(filepre, "_", version, "_", icov, "_", idat, "_", varname,"_", datatime,"_monthly",filesuf, ".nc3", sep="")
        file = paste(outputsdir,trimopt,"/",fname, sep="")
        print(file)
        filepdf = paste(pathplot, fname, ".pdf", sep="")
        pdf(file=filepdf, width=11, height=7)
        ncid = open.nc(con=file, write=FALSE)
        timevec = var.get.nc(ncid, 'time')
        #for (m in 1:length(MON)) {
        for (m in 1:length(timevec)) {
           print(m)
           par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
           for (p in covlist) {
              p = trim(p)
              print(p)
              x = var.get.nc(ncid,p)

              #{x[x==0]=NA}
              plot.grid.continuous(mapz=x[,,m], res=res,colors=drywet(40),
                    xlab="", ylab="", zlim=c(0,7))
              plot(coastsCoarse, add=TRUE, lwd=0.25)
              long_name = att.get.nc(ncid, p, attribute="long_name")
              units = att.get.nc(ncid, p, attribute="units")
              mtext(paste(long_name, " (",units,")", sep=""), cex=0.6)
              #title(paste(p, " (",units,")", sep=""))
           }
              mtext(outer=TRUE, fname, line=1)
              mtext(outer=TRUE, MON[timevec[m]], line=-1)
        }
        dev.off()

        #Plot again with NA for zero's.
        filepdf = paste(pathplot, fname,"_NA", ".pdf", sep="")
        pdf(file=filepdf, width=11, height=7)
        for (m in 1:length(timevec)) {
           par(mfrow=c(4,4), omi=c(0,0.0,.5,0.5), mar=c(1,1,2,2)+0.1)
           for (p in covlist) {
              p = trim(p)
              print(p)
              x = var.get.nc(ncid,p)[,,m]
              x[x==0]=NA
              plot.grid.continuous(mapz=x[,], res=res,colors=drywet(40),
                    xlab="", ylab="", zlim=c(0,7))
              plot(coastsCoarse, add=TRUE, lwd=0.25)
              long_name = att.get.nc(ncid, p, attribute="long_name")
              units = att.get.nc(ncid, p, attribute="units")
              mtext(paste(long_name, " (",units,")", sep=""), cex=0.6)
              #title(paste(p, " (",units,")", sep=""))
           }
              mtext(outer=TRUE, fname, line=1)
              mtext(outer=TRUE, MON[timevec[m]], line=-1)
        }
        dev.off()
        close.nc(ncid)

    }
}

