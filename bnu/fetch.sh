rsync -av gibbs:/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent . --filter='+ *_qxq.nc' --filter='- *.nc' --filter='+ raw/*'


#rsync --dry-run -av gibbs:/home2/rpfische/git/entgvsd1/bnu/lc_lai_ent . --include='*_qxq.nc' --include='trimmed_scaled/*' --include='trimmed/*' --include='maxcrops/*' --include='nocrops/*' --include='purelr/*' --include='raw/*' --exclude='*'

