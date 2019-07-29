# Fetch raw/ and pure/ --- but only _qxq lo-res versions
rsync  -av gibbs:/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent . --filter='+ *_qxq.nc' --filter='- *.nc'

# Fetch everything else
rsync  -av --filter='- /lc_lai_ent/raw' --filter='- /lc_lai_ent/pure' gibbs:/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent . 
