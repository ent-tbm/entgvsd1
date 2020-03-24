#B21b_plots_to_png.R
source("utils.R")
library("imagemagick")

pdflist = read.table("../outputs/plots/filelist.txt")

file = pdflist[1]
im.convert(file, output=paste("../outputs/png/",file, ".png", sep=""), extra.opts="-density 150")

#library(pdftools); pdf_convert(pdflist[1], format = "png", pages = NULL, filenames = NULL, dpi = 300, opw = "", upw = "", verbose = TRUE)


