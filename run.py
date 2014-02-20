#/usr/bin/env python

import os
import sys
import glob
import multiprocessing
import simplejson as json
import yaml
from util import find_files, mkdirs, get_bcd_subset
from bcd_list import get_bcd_list
from bcd_sources import map_bcd_sources
from bcd_phot import get_bcd_phot 
from bcd_phot import write_mean_groups
from bcd_phot import save_single_channel
from bcd_phot import apply_array_location_correction
from bcd_phot import array_location_setup


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

# this is the master project directory containing the data, input file,
# and subdirectory containing the RA/Dec source list files
data_dir = params['data_dir'].rstrip('/')
print('data directory: '+data_dir)

# this is the name of the output dir for pipeline output and temp files
out_dir = params['out_dir'].rstrip('/')
mkdirs(out_dir)
print('output directory: '+out_dir)

# create a dictionary with BCD filenames as the keys, full paths as values
print('creating setup files in output directory...')
bcd_paths = [i for i in find_files(data_dir,bcd_suffix)]
bcd_dict = {i.split('/')[-1]:i for i in bcd_paths}
with open(out_dir+'/bcd_dict.json','w') as w:
	json.dump(bcd_dict,w,indent=' '*4)
# do the same for UNC files
unc_paths = [i for i in find_files(data_dir,unc_suffix)]
unc_dict = {i.split('/')[-1]:i for i in unc_paths}
with open(out_dir+'/unc_dict.json','w') as w:
	json.dump(unc_dict,w,indent=' '*4)

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
			bcd_paths_subset = get_bcd_subset(bcd_dict,aors,ch,hdr)
			unc_paths_subset = get_bcd_subset(unc_dict,aors,ch,hdr)
		else:
			bcd_paths_subset = get_bcd_subset(bcd_dict,aors,ch)
			unc_paths_subset = get_bcd_subset(unc_dict,aors,ch)

		bcd_dict_subset = {i.split('/')[-1]:i for i in bcd_paths_subset}
		unc_dict_subset = {i.split('/')[-1]:i for i in unc_paths_subset}

		bcd_dict_path = work_dir+'/bcd_dict.json'
		unc_dict_path = work_dir+'/unc_dict.json'

		with open(bcd_dict_path,'w') as w:
			json.dump(bcd_dict_subset,w,indent=' '*4)
		print('created: '+bcd_dict_path)

		with open(unc_dict_path,'w') as w:
			json.dump(unc_dict_subset,w,indent=' '*4)
		print('created: '+unc_dict_path)

		metadata = {'name':name, 'data_dir':data_dir, 'work_dir':work_dir,
			'out_dir':out_dir, 'radecfile':radecfile,
			'bcd_dict_path':bcd_dict_path, 'unc_dict_path':unc_dict_path,
			'aors':aors, 'channel':ch, 
			'cbcd':params['cbcd'], 'mask':params['mask'],
			'idl_path': idl_path, 'max_cov': max_cov}
		if is_hdr:
			metadata['hdr'] = hdr
		metadata_path = work_dir+'/metadata.json'
		with open(metadata_path,'w') as w:
			json.dump(metadata,w,indent=' '*4)
		print('created: '+metadata_path)

# loop through the region/channel/hdr file structure just created and 
# read the metadata for each, then call get_bcd_list(meta)
print('associating input sources with BCDs for:')
for work_dir in work_dirs:
	print(work_dir)
	metadata = json.load(open(work_dir+'/metadata.json'))
	get_bcd_list(metadata)

# now run map_to_sources to reverse the mapping such that each BCD has a 
# list of its associated sources
# filepaths = glob.glob(out_dir+'/*/*/*/bcd_list.json')
print('mapping BCDs to their sources...')
filepaths = [i for i in find_files(out_dir,'bcd_list.json')]
pool.map(map_bcd_sources,filepaths)

# now run get_bcd_phot to compute photometry on all the sources 
# source_list_paths = glob.glob(out_dir+'/*/*/*/source_list.json')
print('getting photometry from IDL...')
filepaths = [i for i in find_files(out_dir,'source_list.json')]
pool.map(get_bcd_phot,filepaths)

# now run write_mean_groups to produce ch1/ch2,long/short catalogs
# phot_group_paths = glob.glob('bcd_dirs/*/*/phot_groups.json')
print('combining multiple measurements of sources...')
filepaths = [i for i in find_files(out_dir,'phot_groups.json')]
pool.map(write_mean_groups,filepaths)

# now run save_single_channel to get individual channel/exposure catalogs
print('writing single channel catalogs...')
filepaths = [i for i in find_files(out_dir,'phot_groups_mean.json')]
pool.map(save_single_channel,filepaths)

# now apply the array location dependent correction
print('applying array location correction...')
filepaths = [i for i in find_files(out_dir,'phot_groups_mean.json')]
args_list = array_location_setup(filepaths)
pool.map(apply_array_location_correction,args_list)