#B21b_plots_to_png.R
#Convert pdf to png.

source("utils.R")
library("imagemagick")

args = commandArgs(trailingOnly=TRUE)
print(args)
numargs = length(args)
if (numargs != 1) {
  print ('Usage:  Rscript B21b_plot_to_png.R <path>')
  print ('path = directory containing all pdf files to convert to png.'
  quit()
} 

pdflist = list.files(path=args[1]) #read.table("../outputs/plots/filelist.txt")
print(pdflist)

file = pdflist[1]
im.convert(file, output=paste("../outputs/png/",file, ".png", sep=""), extra.opts="-density 150")

#library(pdftools); pdf_convert(pdflist[1], format = "png", pages = NULL, filenames = NULL, dpi = 300, opw = "", upw = "", verbose = TRUE)


