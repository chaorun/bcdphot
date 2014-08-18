#/usr/bin/env python

import os
import sys
import glob
import multiprocessing
import simplejson as json
import yaml
from util import find_files, mkdirs, get_bcd_subset
from bcd_list import get_bcd_list
from bcd_list import map_bcd_sources
from bcd_phot import get_bcd_phot 
from bcd_phot import apply_array_location_correction
from bcd_phot import calculate_full_uncertainties
from bcd_phot import cull_bad_measurements
from bcd_phot import combine_hdr_catalogs
from bcd_phot import uncorrect_red_sources

# instantiate the pool for parallel processing
ncpus = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=ncpus)
print "using %i CPUs" % ncpus

# read in the YAML setup file containing the list of dicts, where
# each dict corresponds to a region and contains the AOR numbers
# and filepaths containing the list of RA/Dec for that region, as well
# as some parameters such as whether or not the data are HDR mode, and
# whether or not to use the *_cbcd.fits files or the *_bcd.fits files.
infile = sys.argv[1]
setup = yaml.load(open(infile))
print('reading input file: '+infile)

regions = setup['regions']
params = setup['params']

# check to see whether to use CBCD or BCD files 
# (and consequently CBUNC or BUNC files)
if params['cbcd']:
	bcd_suffix = '*_cbcd.fits'
	unc_suffix = '*_cbunc.fits'
else:
	bcd_suffix = '*_bcd.fits'
	unc_suffix = '*_bunc.fits'

# check to see whether the data are HDR mode or not. if they are then
# there will be additional directory structure in the output directory.
is_hdr = params['hdr']

# get the path to the IDL executable from the params dict
idl_path = params['idl_path']

# get the number corresponding to the maximum coverage of any source in
# the input source lists. an unnecessarily large value will make the 
# pipeline slower than it has to be, so this value should be set to the
# appropriate value depending on the observing strategy. if the value is
# too low then the pipeline will not take advantage of additional data
# available for some sources.
max_cov = params['max_cov']

# get the minimum SNR value for an individual measurement to be kept
min_snr = params['min_snr']

# get the maximum distance from input RA/Dec to centroid position [arcsec]
max_dist = params['max_dist']

# this is the master project directory containing the data, input file,
# and subdirectory containing the RA/Dec source list files
data_dir = params['data_dir'].rstrip('/')
print('data directory: '+data_dir)

# this is the name of the output dir for pipeline output and temp files
out_dir = params['out_dir'].rstrip('/')
mkdirs(out_dir)
print('output directory: '+out_dir)

# if out_dir not empty, assume bcd_dict.json and unc_dict.json exist 
# and are valid (this can save time if running the same region
# more than once, i.e. during testing)
if len(os.listdir(out_dir)) is 0:

	# create a dictionary with BCD filenames as the keys, full paths as values
	print('creating setup files in output directory...')
	bcd_paths = list(find_files(data_dir, bcd_suffix))
	bcd_dict = {i.split('/')[-1]:i for i in bcd_paths}
	with open(out_dir+'/bcd_dict.json','w') as w:
		json.dump(bcd_dict, w, indent=' '*4)
	# do the same for UNC files
	unc_paths = list(find_files(data_dir, unc_suffix))
	unc_dict = {i.split('/')[-1]:i for i in unc_paths}
	with open(out_dir+'/unc_dict.json','w') as w:
		json.dump(unc_dict, w, indent=' '*4)
	# do the same for MSK files
	msk_paths = list(find_files(data_dir, '*_bimsk.fits'))
	msk_dict = {i.split('/')[-1]:i for i in msk_paths}
	with open(out_dir+'/msk_dict.json','w') as w:
		json.dump(msk_dict, w, indent=' '*4)

else:

	# read in existing bcd_dict.json, unc_dict.json
	bcd_dict = json.load(open(out_dir+'/bcd_dict.json'))
	unc_dict = json.load(open(out_dir+'/unc_dict.json'))
	msk_dict = json.load(open(out_dir+'/msk_dict.json'))

# loop through the source lists (radecfiles) and create output directory
# structure and metadata files used throughout the rest of the pipeline
print('creating directory structure and metadata in output directory...')
work_dirs = []
for region in regions:

	name = region['name']
	radecfiles = region['radec']
	aors = region['aors']

	for radecfile in radecfiles:
		f = radecfile.split('/')[-1]
		if is_hdr:
			ch, hdr = f.split('_')[1:3]
			work_dir = '/'.join([out_dir,name,ch,hdr])
		else:
			ch = f.split('_')[1]
			work_dir = '/'.join([out_dir,name,ch])
		mkdirs(work_dir)
		work_dirs.append(work_dir)

		# make subset bcd_dict/unc_dict in each work_dir
		if is_hdr:
			bcd_paths_subset = get_bcd_subset(bcd_dict, aors, ch, hdr)
			unc_paths_subset = get_bcd_subset(unc_dict, aors, ch, hdr)
			msk_paths_subset = get_bcd_subset(msk_dict, aors, ch, hdr)
		else:
			bcd_paths_subset = get_bcd_subset(bcd_dict, aors, ch)
			unc_paths_subset = get_bcd_subset(unc_dict, aors, ch)
			msk_paths_subset = get_bcd_subset(msk_dict, aors, ch)

		bcd_dict_subset = {i.split('/')[-1]:i for i in bcd_paths_subset}
		unc_dict_subset = {i.split('/')[-1]:i for i in unc_paths_subset}
		msk_dict_subset = {i.split('/')[-1]:i for i in msk_paths_subset}

		bcd_dict_path = work_dir+'/bcd_dict.json'
		unc_dict_path = work_dir+'/unc_dict.json'
		msk_dict_path = work_dir+'/msk_dict.json'

		with open(bcd_dict_path,'w') as w:
			json.dump(bcd_dict_subset, w, indent=' '*4)
		print('created: '+bcd_dict_path)

		with open(unc_dict_path,'w') as w:
			json.dump(unc_dict_subset, w, indent=' '*4)
		print('created: '+unc_dict_path)

		with open(msk_dict_path,'w') as w:
			json.dump(msk_dict_subset, w, indent=' '*4)
		print('created: '+msk_dict_path)

		metadata = {'name':name, 'data_dir':data_dir, 'work_dir':work_dir,
			'out_dir':out_dir, 'radecfile':radecfile,
			'bcd_dict_path':bcd_dict_path,
			'unc_dict_path':unc_dict_path,
			'msk_dict_path':msk_dict_path,
			'aors':aors, 'channel':ch, 
			'cbcd':params['cbcd'], 'mask':params['mask'],
			'idl_path': idl_path, 'max_cov': max_cov,
			'min_snr':min_snr, 'max_dist':max_dist,
			'sigma_clip':params['sigma_clip'],
			'centroid':params['centroid']
			}
		if is_hdr:
			if ch == '1':
				metadata['long_cutoff'] = params['long_cutoff_ch1']
				metadata['short_cutoff'] = params['short_cutoff_ch1']
			elif ch == '2':
				metadata['long_cutoff'] = params['long_cutoff_ch2']
				metadata['short_cutoff'] = params['short_cutoff_ch2']
			metadata['hdr'] = hdr
		metadata_path = work_dir+'/metadata.json'
		with open(metadata_path,'w') as w:
			json.dump(metadata, w, indent=' '*4)
		print('created: '+metadata_path)

# if bcd_list.json files exist, assume valid and don't run get_bcd_list
if len(list(find_files(out_dir, 'bcd_list.json'))) is 0:

	# loop through the region/channel/hdr file structure just created and 
	# read the metadata for each, then call get_bcd_list(meta)
	print('associating input sources with BCDs for:')
	for work_dir in work_dirs:
		print(work_dir)
		metadata = json.load(open(work_dir+'/metadata.json'))
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
# print('culling bad measurements...')
# filepaths = list(find_files(out_dir, 'phot_groups.json'))
# pool.map(cull_bad_measurements, filepaths)

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
	# insert non-hdr sigma clip code here
	pass
