#EntGVSD_efnk.R
source("utils.R")

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
    path = "/Users/nkiang/NancyResearch/GISS/Models/entgvsd1/entgvsd1_1/bnu/"
    #Choose date of output
    #dateout = "2019-07-19"
    #dateout = "2019-08-22"
    #dateout = "2019-08-26"
    #dateout = "2019-08-29"
    dateout = "2019-09-03"
    #entlclaidir = paste(dateout,"/lc_lai_ent.v7/", sep="")
    #entlclaidir = paste(path, dateout,"/lc_lai_ent/", sep="")
    entlclaidir = paste(path,"/lc_lai_ent/", sep="")
        #If in R gui:
    #pathplot=paste(entlclaidir, "/PLOTS/", sep="")
        #If in git clone:
    pathplot="PDFPLOTS/"
        dir.create(pathplot)
}
#For Nancy on gibbs:
if (FALSE) {
    entlclaidir = "/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent/"
    pathplot = "PLOTS/" 
    # Create output directory if it doesn't already exist
    dir.create(pathplot)
}
    
    
#For running within the BNU/ directory as part of the generation process:
if (TRUE) {
    entlclaidir = "lc_lai_ent/"
    pathplot=paste(entlclaidir, "plots/", sep="")
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
  
    map.entgvsd.steps(entlclaidir, res=res, enttyp,varname="lc", trimopt=trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="laimax", trimopt=trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(entlclaidir, res=res,  enttyp, varname="hgt", trimopt=trimopt, filepre, datatime, version,  icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    for (d in c( "2004_017", "2004_201")) {
        map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version, icov, idat, filesuf=filesuf,do.pdf=do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }
    for (d in paste(datatime, "_", MON, sep="")) {
        map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version, icov, idat, filesuf=filesuf,do.pdf=do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }

	#Layered check (only for ent17):
    map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="laimax_err", trimopt=trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = FALSE) #This is by PFT

    #Single layer checks
    # Monthly
    for (d in paste(datatime, "_", MON, sep="")) {
    	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    	
	}
	
	# Totals
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_dompft", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    	
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_dompftlc", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	
 	
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.

   #map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lclaimax_err", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
    
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.
 
 	#map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lchgt_err", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
  	
  	#Checksum diff maps
  	res="HXH"
  	checksumdir = paste(entlclaidir, "checksum/", sep="")
  	checksuff = paste("_diff", filesuf, sep="")
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	for (d in paste(datatime, "_", MON, sep="")) { 	 
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
    for (d in c( "2004_017", "2004_201")) {
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 

	
	
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
    
    map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="lc",  trimopt=trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="laimax",  trimopt=trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    map.entgvsd.steps(entlclaidir, res=res,  enttyp=enttyp, varname="hgt",  trimopt=trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    #for (d in c( "2004_Jan", "2004_Jul")) {
    for (d in paste(datatime, "_", MON, sep="")) {
        map.entgvsd.steps(entlclaidir, res=res, enttyp=enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version=version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }
 
    for (d in paste(datatime, "_", MON, sep="")) {
        map.entgvsd.steps(entlclaidir, res=res, enttyp, varname="lai", trimopt=trimopt, filepre, datatime=d, version, icov, idat, filesuf=filesuf,do.pdf=do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    }

    #Single layer checks
    # Monthly
    for (d in paste(datatime, "_", MON, sep="")) {
    	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    
	}

		
	# Totals
	# Dominant PFT 
	fname = paste(trimopt,"/", filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	domlc = Ent_calc_domlc(file=paste(entlclaidir, fname, sep=""), enttyp)
	fnameout = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc_dompft", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
	pdf(file=paste(pathplot, fnameout, ".pdf", sep=""), width=8, height=5)
	Ent_dompft_plot(lctype=domlc, numpft=16, res=res, legend.cex=0.6, Entcolors=Entcolors16, if.new=FALSE)
   	dev.off()
   	
	if (FALSE) {
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_dompft", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    	
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_dompftlc", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 	} #if (FALSE)
 	
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.

	
  	#Checksum diff maps
  	res = "HXH"
  	checksumdir = paste(entlclaidir, "checksum/",  sep="")
  	checksuff = paste("_diff", filesuf, sep="")
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	for (d in paste(datatime, "_", MON, sep="")) { 	 
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
    for (d in c( "2004_017", "2004_201")) {
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 

} #pure
 
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
	        map.entgvsd.steps(entlclaidir, res=res, varname, enttyp=enttyp, trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
	        #map.entgvsd.steps(entlclaidir, res=res, varname="lc", enttyp=enttyp, trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    	    #map.entgvsd.steps(entlclaidir, res=res, varname="laimax", enttyp=enttyp, trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
        	#map.entgvsd.steps(entlclaidir, res=res, varname="hgt",enttyp=enttyp, trimopt, filepre, datatime, version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    	}
    	
    	for (d in paste(datatime, "_", MON, sep="")) {
    	    map.entgvsd.steps(entlclaidir, res=res, enttyp=enttyp, varname="lai", trimopt, filepre, datatime=d, version=version, icov, idat, filesuf=filesuf,do.pdf = do.pdf, pathplot=pathplot, do.checksum = do.checksum)
    	}
    
   		#Single layer checks
    	# Monthly
    	for (d in paste(datatime, "_", MON, sep="")) {
    		map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    
		}
	    	
		# Totals
		fname = paste(trimopt,"/", filepre, "_",version, "_",icov, "_",idat, "_","lc", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
		domlc = Ent_calc_domlc(file=paste(entlclaidir, fname, sep=""), enttyp)
		fnameout = paste(filepre, "_",version, "_",icov, "_",idat, "_","lc_dompft", "_",datatime,"_ann",  "_", trimopt, filesuf, ".nc", sep="")
		pdf(file=paste(pathplot, fnameout, ".pdf", sep=""), width=8, height=5)
		Ent_dompft_plot(lctype=domlc, numpft=16, res=res, legend.cex=0.6, Entcolors=Entcolors16, if.new=FALSE)
   		dev.off()

		if (FALSE) {
   	 	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_dompft", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)    
    	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_npftgrid", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 		map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lc_dompftlc", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 		} #if (FALSE)
 		
    	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with laimax.
    	map.entgvsd.check.misc(entlclaidir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot)  #This is also plotted with hgt.
 	
 	#Checksum diff maps
  	checksumdir = paste(entlclaidir, "checksum/", sep="")
  	checksuff = "_diff"
  	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclaimax_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	for (d in paste(datatime, "_", MON, sep="")) { 	 
  	 	map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lclai_checksum", trimopt, filepre, datatime=d,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
	}
 	 map.entgvsd.check.misc(checksumdir, res, enttyp=enttyp, varnamecheck="lchgt_checksum", trimopt, filepre, datatime,  version, icov, idat, filesuf=checksuff, add.new=FALSE, do.pdf = do.pdf, pathplot=pathplot) 
 		
   } #trimopts
} #trimming
