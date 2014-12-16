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

setup_output_dirs(setup)

# if bcd_list.json files exist, assume valid and don't run get_bcd_list
if len(list(find_files(out_dir, 'bcd_list.json'))) is 0:

	# loop through the region/channel/hdr file structure just created and 
	# read the metadata for each, then call get_bcd_list(meta)
	print('associating input sources with BCDs for:')
	for metafile in find_files(out_dir, 'metadata.json'):
		print(metafile)
		metadata = json.load(open(metafile))
		get_bcd_list(metadata)

# if source_list.json files exist, assume valid and don't run map_bcd_sources
if len(list(find_files(out_dir, 'source_list.json'))) is 0:
	# now run map_to_sources to reverse the mapping such that each BCD has a 
	# list of its associated sources
	# filepaths = glob.glob(out_dir+'/*/*/*/bcd_list.json')
	print('mapping BCDs to their sources...')
	filepaths = list(find_files(out_dir, 'bcd_list.json'))
	pool.map(map_bcd_sources, filepaths)

# run get_bcd_phot to compute photometry on all the sources 
# source_list_paths = glob.glob(out_dir+'/*/*/*/source_list.json')
print('getting photometry from IDL...')
filepaths = list(find_files(out_dir, 'source_list.json'))
pool.map(get_bcd_phot, filepaths)


# cull bad individual measurements from the groups of measurements
# of each source using SNR and proximity cutoffs
if setup['params']['snr_dist_cull']:
	print('culling bad measurements...')
	filepaths = list(find_files(out_dir, 'phot_groups.json'))
	pool.map(cull_bad_measurements, filepaths)

# apply array location correction to phot_groups.json files
print('applying array location correction...')
# filepaths = list(find_files(out_dir, 'phot_groups_culled.json'))
filepaths = list(find_files(out_dir, 'phot_groups.json'))
pool.map(apply_array_location_correction, filepaths)

# un-correct photometry of red sources, by matching ch1/ch2
if is_hdr:
	filepaths_long = filter(lambda x: 'long' in x, 
		find_files(out_dir, 'phot_groups_arrayloc.json'))
	filepaths_short = filter(lambda x: 'short' in x, 
		find_files(out_dir, 'phot_groups_arrayloc.json'))
	ch1long = filter(lambda x: '/1/' in x, filepaths_long)
	ch2long = filter(lambda x: '/2/' in x, filepaths_long)
	ch1short = filter(lambda x: '/1/' in x, filepaths_short)
	ch2short = filter(lambda x: '/2/' in x, filepaths_short)
	pool.map(uncorrect_red_sources, zip(ch1long, ch2long))
	pool.map(uncorrect_red_sources, zip(ch1short, ch2short))
else:
	filepaths = list(find_files(out_dir, 'phot_groups_arrayloc.json'))
	ch1 = filter(lambda x: '/1/' in x, filepaths)
	ch2 = filter(lambda x: '/2/' in x, filepaths)
	pool.map(uncorrect_red_sources, zip(ch1, ch2))

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