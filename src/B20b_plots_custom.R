# Plots single file results of EntGVSD, by specifying file type and name at the command line.
# 
# Author: Nancy Kiang
# 
# To run:
#   Rscript B20b_plots_custome.R NPFT filetype path filename
# where:
#   [NPFT] = 16 | 17
#   [filetype] = lc | laimax | laimonth | laidoy | hgt | lc_checksum | lclai_checksum | lchgt_checksum

source("utils.R")

#******************************************************************************************
# Step 1:  EDIT THIS SECTION FOR YOUR PATHS AND TO SET UP A plots DIRECTORY.

#Set paths:

#Local directory customization, example.
if (FALSE) {
    outputsdir = "../outputs/"
    pathplot = "../outputs/plots/" 
    # Create output directory if it doesn't already exist
    dir.create("../outputs/")
    dir.create(pathplot)
}
    
#DEFAULT: For running within the src/ directory as part of the generation process:
if (TRUE) {
    pathplot="../myplots"
    # Create output directory if it doesn't already exist
    dir.create(pathplot)
}

#*****************************************************************************************
B20b_usage = function() {
  print("Usage: Rscript B20b_plots_custom.R NPFT filetype path filename <vname>")
  print("<vname> required if filetype is soilalbedo.")
  print("npft = 16 | 17")
  print("filetype = lc | laimax | laimonth | laidoy | hgt | lc_checksum | lclai_checksum | lchgt_checksum")
  print(c("lc_checksum", "lclai_checksum", "lclaimax_checksum", "lchgt_checksum") )
  print(c("lc_dompft", "lc_nftpgrid", "lc_domlc"))
  print("soilalbedo")
}

#*****************************************************************************************
args = commandArgs(trailingOnly=TRUE)
print(args)
numargs = length(args)
if (numargs < 4) {
  B20b_usage()
  quit()
} 
npft = as.numeric(args[1])
ftype = args[2]
path = args[3]
fname = args[4]

if (npft==17) { #ent17
   ncov = npft + 3
} else if (npft==16) { #pure onwards
   ncov = npft + 2
}

file = paste(path, "/", fname, sep="")
print(file)

if (ftype=="soilalbedo" & numargs<5) {
   B20b_usage()
   quit()
}

if (numargs==5) {
     vname = args[5]
}

ftype3D = c("lc", "laimax", "laimonth", "laidoy", "hgt")
ftypecs = c("lc_checksum", "lclai_checksum", "lclaimax_checksum", "lchgt_checksum") 
ftypeother = c("lc_dompft", "lc_nftpgrid", "lc_domlc")
ftypesoil = c("soilalbedo")
filetypes = c(ftype3D, ftypecs, ftypeother, ftypesoil)

if (is.na(match(ftype, filetypes))) {
   print(c("No such file type: ", ftype))
   B20b_usage()
   quit()
}

#******************************************************************************************
# Step 2: Get resolution
nci = open.nc(con=file)
lon = dim.inq.nc(nci, 'lon')$length
res = res.from.IM(IM=lon)
print(c(lon, res))

#******************************************************************************************
# Step 3:  Select output to screen or to pdf.  Select whether to output checksums.
do.pdf = TRUE

   if (do.pdf) {
        filepdf = paste(pathplot,"/",fname, ".pdf", sep="")
        pdf(file=filepdf, width=11, height=7)
    } else if (add.new) {
        quartz(width=11,height=7) #Open plotting window
    } else {
        print("Printing to screen.")
    }

#******************************************************************************************
# Step 4:  Do plots.
#

if (!is.na(match(ftype, ftype3D))) {
        par(mfrow=c(4,5), omi=c(0,0.0,1.0,0.5), mar=c(1,1,2,2)+0.1)
	#lcnum=c(2,4,6,7:15)
	lcnum=c(1:ncov)
        map.EntGVSD.v1.1(filelc=NULL, file=file, res=res, varpre="", varname=ftype, lcnum=lcnum, colors=giss.palette.nowhite(40), type="any", zlim=c(0,1), if.zeroNA=TRUE, titletype=1)

}

if (!is.na(match(ftype, ftypecs))) {

	#lcnum=c(2,4,6,7:15)
	lcnum=c(1:ncov)
        map.EntGVSD.v1.1(filelc=NULL, file=file, res=res, varpre="", varname=ftype, lcnum=lcnum, colors=giss.palette.nowhite(40), type="any", zlim=c(0,1), if.zeroNA=TRUE, titletype=1)

} 


if (!is.na(match(ftype, ftypecs))) {
        varnamecheck = ftype
	trimopt = npft

        if (varnamecheck ==  "lc_dompft" | varnamecheck == "lc_domlc") {
                zlim = c(min(enttyp), max(enttyp))
                if (trimopt=="ent17") {
                        color = Entrgbhex(Entcolors17b[1:20,])
                        leg=Entcolors17b[1:20, "lc_type"]
                } else {
                        color = Entrgbhex(Entcolors16)
                        leg=Entcolors16[,"lc_type"]
                }
                restime = "_ann"
                if.cat=TRUE

        } else if (varnamecheck == "lc_npftgrid" ) {
                if (FALSE) { #categorical colors do not work
                  zlim = c(0,17)
                  color = giss.palette.nowhite(18)
                  color[1] = rgb(0.5, 0.5, 0.5)
                  leg=zlim[1]:zlim[2]
                  restime="_ann"
                  if.cat=TRUE
                }
                zlim = c(0,17)
                color = giss.palette.nowhite(18)
                color[1] = rgb(0.5, 0.5, 0.5)
		leg=zlim[1]:zlim[2]
                restime="_ann"
                if.cat=FALSE
        } else if (varnamecheck == "lc_dompftlc" | varnamecheck == "lc_checksum") {
                zlim = c(0,1)
 	        color = giss.palette.nowhite(11)
                color[1] = rgb(0.5, 0.5, 0.5)
                leg=(0:10)/10
                restime = "_ann"
		if.cat=FALSE
	}
        map.GCM (file, varname=varnamecheck, res=res,colors=color,  zlim=zlim, if.cat=if.cat, if.zeroNA=FALSE, titletype=1)
}

if (ftype=="soilalbedo" & numargs==5) {
     zlim = c(0,1)
     color = giss.palette.nowhite(40)

     mapz = var.get.nc(nci, vname)	
     plot.grid.continuous(mapz=mapz, res=res, colors=color, legend.lab=NULL, xlab="longitude", ylab="latitude", titletext="", zlim=zlim, ADD=FALSE, if.fill=TRUE, if.coasts=TRUE, ask=TRUE)
       
}

if (do.pdf) {
   dev.off()
}


