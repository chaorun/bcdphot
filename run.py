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
from bcd_phot import write_mean_groups
from bcd_phot import save_single_channel
from bcd_phot import apply_array_location_correction


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
	proj_dir = sys.argv[2].rstrip('/')
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
			# make subset bcd_dict/unc_dict in each work_dir
			bcd_paths_subset = get_bcd_subset(bcd_dict,aors,ch,hdr)
			unc_paths_subset = get_bcd_subset(unc_dict,aors,ch,hdr)
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
			metadata = {'name':name, 'proj_dir':proj_dir, 'work_dir':work_dir,
				'out_dir':out_dir, 'radecfile':radecfile,
				'bcd_dict_path':bcd_dict_path, 'unc_dict_path':unc_dict_path,
				'aors':aors, 'channel':ch, 'hdr':hdr}
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
	ch1, ch2, out_paths = [], [], []
	for i in range(len(filepaths)):
		work_dir = filepaths[i].split('/phot_groups_mean.json')[0]
		metadata = json.load(open(work_dir+'/metadata.json'))
		channel = metadata['channel']
		if channel == '1':
			ch1.append(filepaths[i])
		elif channel == '2':
			ch2.append(filepaths[i])
	# make sure the ith elements of ch1 and ch2 are for the same region/exptime
	ch1.sort()
	ch2.sort()
	# construct file paths for matched ch1/ch2 catalogs
	for i in range(len(ch1)):
		work_dir = ch1[i].split('/phot_groups_mean.json')[0]
		metadata = json.load(open(work_dir+'/metadata.json'))
		out_dir = '/'.join( [ metadata['out_dir'], metadata['name'] ] )
		spl1 = ch1[i].split('bcdphot_out/')[1].split('/')
		spl2 = ch2[i].split('bcdphot_out/')[1].split('/')
		if spl1[0] != spl2[0] or spl1[2] != spl2[2]:
			print('error: ch1/ch2 catalogs don\'t match')
		else:
			out_name = '_'.join([spl1[0],spl1[2],'matched_catalog.txt'])
			out_paths.append('/'.join([out_dir,out_name]))
	args_list = zip(ch1, ch2, out_paths)
	pool.map(apply_array_location_correction,args_list)