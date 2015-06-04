#/usr/bin/env python

import os
import sys
import multiprocessing
import simplejson as json
import yaml
from util import find_files, setup_output_dirs
from bcd_list import get_bcd_list
from bcd_list import map_bcd_sources
from bcd_phot import get_bcd_phot 
from bcd_phot import apply_array_location_correction
from bcd_phot import calculate_full_uncertainties
from bcd_phot import cull_bad_measurements
from bcd_phot import combine_hdr_catalogs
from bcd_phot import uncorrect_red_sources
from bcd_phot import sigma_clip_non_hdr
from bcd_phot import make_2ch_catalogs

# instantiate the pool for parallel processing
ncpus = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=ncpus)
print "using %i CPUs" % ncpus

infile = sys.argv[1]
setup = yaml.load(open(infile))
print('reading input file: '+infile)

is_hdr = setup['params']['hdr']
out_dir = setup['params']['out_dir'].rstrip('/')

# calculate the full uncertainties from both systematic
# and photometric processes
print('calculating full uncertainties...')
filepaths = list(find_files(out_dir, 'phot_groups_arrayloc_cor.json'))
pool.map(calculate_full_uncertainties, filepaths)

# create catalogs - since this step operates on fully corrected data
# with full uncertainties, this is where we apply a basic sigma-clip
print('creating catalogs...')
if is_hdr:
	filepaths_long = list(find_files(out_dir, '*long_catalog.txt'))
	filepaths_short = list(find_files(out_dir, '*short_catalog.txt'))
	ch1long = [i for i in filepaths_long if '_1_' in i]
	ch2long = [i for i in filepaths_long if '_2_' in i]
	ch1short = [i for i in filepaths_short if '_1_' in i]
	ch2short = [i for i in filepaths_short if '_2_' in i]
	pool.map(combine_hdr_catalogs, zip(ch1long, ch1short))
	pool.map(combine_hdr_catalogs, zip(ch2long, ch2short))
else:
	filepaths = list(find_files(out_dir, '*catalog.txt'))
	pool.map(sigma_clip_non_hdr, filepaths)
	filepaths_ch1 = list(find_files(out_dir, '*_1_catalog_sigclip.txt'))
	filepaths_ch2 = list(find_files(out_dir, '*_2_catalog_sigclip.txt'))
	pool.map(make_2ch_catalogs, zip(filepaths_ch1, filepaths_ch2))

pool.close()
