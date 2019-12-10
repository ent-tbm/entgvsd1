The original 5min files were transformed to NetCDF with:

giss2nc --input-file Monfreda_crops_5min.bin --output-file Monfreda_crops_5min.nc --endian=big --names='c3crop,c4crop,herbcrop,shrubcrop,treecrop,c3c4tot,herbshrubtreetot'

giss2nc --endian=big --input-file=Monfreda_crops_5min_norm.bin --output-file=Monfreda_crops_norm_5min.nc --names='c3crop,c3crop2,c3c4tot,herbcrop,shrubcrop,treecrop,herbshrubtreetot,c3cropmh,c4cropmh,c3herbcropmh,c4herbcropmh,shrubtreetot,c3c4tot'

giss2nc --endian=big --input-file=Monfreda_crops_5min_normA.bin --output-file=Monfreda_crops_normA_5min.nc --names='c3crop,c4crop,c4cropmh,shrupcrop,treecrop,c3crop2,c4crop,c3c4tot,herbcrop,shrubcrop,treecrop,herbshrubtreecrop'

See here for giss2nc program, which converts "old-format" GISS Fortran binary files to NetCDF format.

https://github.com/citibeth/icebin/blob/develop/modele/giss2nc.cpp
