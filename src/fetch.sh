# Sample rsync commands to fetch plots and SMALL output files to local
# computer where they can be viewed.

# Fetch raw/ and pure/ --- but only _qxq lo-res versions
rsync  -av gibbs:/home2/rpfische/git/entgvsd1/outputs . --filter='+ *_forplot.nc' --filter='- *.nc' --filter='- *.tar.gz' --filter='- *.pdf'

# Fetch everything else
rsync  -av --filter='- /lc_lai_ent/ent17' --filter='- /lc_lai_ent/carrer' --filter='- /lc_lai_ent/pure' gibbs:/home2/rpfische/git/entgvsd1/outputs . 

# Alternatively, just fetch plots
rsync -av gibbs:/home2/rpfische/git/entgvsd1/outputs/png lc_lai_ent
