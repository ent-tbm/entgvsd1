# Fetch raw/ and pure/ --- but only _qxq lo-res versions
rsync  -av gibbs:/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent . --filter='+ *_qxq.nc' --filter='+ *_plot.nc' --filter='- *.nc' --filter='- *.tar.gz' --filter='- *.pdf'

## Fetch everything else
rsync  -av --filter='- /lc_lai_ent/ent17' --filter='- /lc_lai_ent/carrer' --filter='- /lc_lai_ent/pure' gibbs:/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent . 

# Fetch plots
#rsync -av gibbs:/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent/png lc_lai_ent
