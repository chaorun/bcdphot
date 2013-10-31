#/usr/bin/env python

import os
import sys
import glob
import multiprocessing
import simplejson as json
from util import find_files, mkdirs, get_bcd_subset
from bcd_list import get_bcd_list
from bcd_sources import map_bcd_sources
from bcd_phot import get_bcd_phot

if __name__ == "__main__":

	# instantiate the pool for parallel processing
	ncpus = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(processes=ncpus)
	print "using %i CPUs" % ncpus

	# this is the JSON input file containing the list of dicts, where
	# each dict corresponds to a region and contains the AOR numbers
	# and filepaths containing the list of RA/Dec for that region
	infile = sys.argv[1]
	regions = json.load(open(infile))
	print('reading input file: '+infile)
	
	# this is the master project directory containing the data, input file,
	# and subdirectory containing the RA/Dec source list files
	proj_dir = sys.argv[2].strip('/')
	print('project directory: '+proj_dir)

	# create output dir for all pipeline output and temporary files
	out_dir = '/'.join([proj_dir,'bcdphot_out'])
	mkdirs(out_dir)
	print('output directory: '+out_dir)

	# create a dictionary with BCD filenames as the keys, full paths as values
	print('creating setup files in project directory...')
	bcd_paths = [i for i in find_files(proj_dir,'*cbcd.fits')]
	bcd_dict = {i.split('/')[-1]:i for i in bcd_paths}
	with open(out_dir+'/bcd_dict.json','w') as w:
		json.dump(bcd_dict,w,indent=' '*4)

	# do the same for UNC files
	unc_paths = [i for i in find_files(proj_dir,'*cbunc.fits')]
	unc_dict = {i.split('/')[-1]:i for i in unc_paths}
	with open(out_dir+'/unc_dict.json','w') as w:
		json.dump(unc_dict,w,indent=' '*4)

	# loop through the source lists (radecfiles) and create output directory
	# structure and metadata files used throughout the rest of the pipeline
	print('creating directory structure and metadata in output directory...')
	work_dirs = []
	for region in regions:
		name = region['name']
		radecfiles = ['/'.join([proj_dir,i]) for i in region['radec']]
		aors = region['aors']
		for radecfile in radecfiles:
			f = radecfile.split('/')[-1]
			ch, hdr = f.split('_')[1:3]
			work_dir = '/'.join([out_dir,name,ch,hdr])
			mkdirs(work_dir)
			work_dirs.append(work_dir)
			bcd_paths = get_bcd_subset(bcd_dict,aors,ch,hdr)
			unc_paths = get_bcd_subset(unc_dict,aors,ch,hdr)
			metadata = {'name':name, 'proj_dir':proj_dir, 'work_dir':work_dir,
				'out_dir':out_dir, 'radecfile':radecfile, 'bcd_paths':bcd_paths,
				'unc_paths':unc_paths, 'aors':aors, 'channel':ch, 'hdr':hdr}
			with open(work_dir+'/metadata.json','w') as w:
				json.dump(metadata,w,indent=' '*4)

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

	# now run get_catalog to produce ch1/ch2,long/short catalogs
	# phot_group_paths = glob.glob('bcd_dirs/*/*/phot_groups.json')
	print('combining multiple measurements of sources...')
	filepaths = [i for i in find_files(out_dir,'phot_groups.json')]
	pool.map(save_catalog,filepaths)

	# now apply the array location dependent correction
	# print('applying array location correction...')
	