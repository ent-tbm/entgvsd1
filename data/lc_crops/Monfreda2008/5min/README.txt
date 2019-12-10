The original data files are from:
 Monfreda, C., N. Ramankutty, and J. A. Foley (2008), Farming the planet: 2. Geographic distribution of crop areas, yields, physiological types, and net primary production in the year 2000, Global Biogeochemical Cycles, 22(1), 19, doi:Gb1022

These are in originally at 5 minutes x 5 minutes spatial resolution, for one year ascii and netcdf format, distinguishing C3 and C4 crops, growth form (tree, shrub, herb), and legumes.  The original data are not in units of strictly cover fraction but count multi-harvest crops more than once during the year, such that total crop cover could sum to great than one.  To derive crop cover fractions for GCM grids, these multi-harvest cover were normalized to 1.  


Any .bin files generated in earlier versions of processing are "GISS layer format" binary files from post-processing with the program hntr4_monfreda2008.f.  So, these are not original data.  Elizabeth states:  The original 5min files were transformed to NetCDF with:

giss2nc --input-file Monfreda_crops_5min.bin --output-file Monfreda_crops_5min.nc --endian=big --names='c3crop,c4crop,herbcrop,shrubcrop,treecrop,c3c4tot,herbshrubtreetot'

giss2nc --endian=big --input-file=Monfreda_crops_5min_norm.bin --output-file=Monfreda_crops_norm_5min.nc --names='c3crop,c3crop2,c3c4tot,herbcrop,shrubcrop,treecrop,herbshrubtreetot,c3cropmh,c4cropmh,c3herbcropmh,c4herbcropmh,shrubtreetot,c3c4tot'

giss2nc --endian=big --input-file=Monfreda_crops_5min_normA.bin --output-file=Monfreda_crops_normA_5min.nc --names='c3crop,c4crop,c4cropmh,shrupcrop,treecrop,c3crop2,c4crop,c3c4tot,herbcrop,shrubcrop,treecrop,herbshrubtreecrop'

See here for giss2nc program, which converts "old-format" GISS Fortran binary files to NetCDF format.

https://github.com/citibeth/icebin/blob/develop/modele/giss2nc.cpp
