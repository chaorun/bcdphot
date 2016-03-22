import os, sys
sys.path.append(os.path.abspath('..'))
from bcdphot.xsc import run_xsc_phot
from bcdphot.hdr import cleanup_sub
from bcdphot.hdr import merge_subarray
from bcdphot.sdss import match_sdss
from bcdphot.sdss import plot_sdss
from bcdphot.wise import match_wise
from bcdphot.wise import plot_wise
from bcdphot.wise import plot_spz_vs_wise
from bcdphot.wise import plot_spz_vs_wise_sdss_class


bcdphot_out_path = '/Users/jlivings/data/spitzer_cirb_data/ace7/bcdphot_out_final'
mosaic_path = '/Users/jlivings/data/spitzer_cirb_data/mopex'
run_xsc_phot(bcdphot_out_path, mosaic_path)

vg_dir = '/Users/jlivings/data/spitzer_cirb_data/subarray/finalphot10'
cleanup_sub(vg_dir)

vg_dir += '_clean'
merge_subarray(vg_dir, bcdphot_out_path)

cat_path = '/Users/jlivings/data/spitzer_cirb_data/subarray/finalphot10_plots_catalogs'
match_sdss(cat_path)
plot_sdss(cat_path)

match_wise(cat_path, sdss=False)
match_wise(cat_path, sdss=True)
plot_wise(cat_path)
plot_spz_vs_wise(cat_path, plot_style='scatter')
plot_spz_vs_wise(cat_path, plot_style='hexbin')
plot_spz_vs_wise(cat_path, plot_style='hist2d')
plot_spz_vs_wise_sdss_class(cat_path, plot_style='scatter')
plot_spz_vs_wise_sdss_class(cat_path, plot_style='hexbin')
plot_spz_vs_wise_sdss_class(cat_path, plot_style='hist2d')
