#/usr/bin/env python

import os
import sys
import glob
import multiprocessing
import simplejson as json
from util 			import find_files
from get_bcd_list	import get_bcd_list
from map_bcd_sources import map_bcd_sources
from get_bcd_phot	import get_bcd_phot

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
	
	# this is the download directory containing all the data downloaded 
	# from the SSC (i.e. dir containing r43432192, r43420416, etc.)
	data_dir = sys.argv[2]

	# decorator for handling exceptions when dir already exists.
	# unnecessary, just for fun
	def try_os(f):
		def wrapper(arg):
			try:
				f(arg)
			except OSError as e:
				pass
		return wrapper

	@try_os
	def mkdir(dir):
		os.mkdir(dir)

	# create output dir
	out_dir = 'bcdphot_out'
	mkdir(out_dir)

	# run get_bcd_list() to associate individual BCDs to the input sources.
	# this function uses multiprocessing by default so just loop through
	# the source lists, and create output directory structure for pipeline
	for region in regions:
		name = region['name']
		# make subdirectory for region
		mkdir('/'.join([out_dir,name]))
		radecfiles = region['radec']
		aors = region['aors']
		for radecfile in radecfiles:
			f = radecfile.split('/')[-1]
			ch, hdr = f.split('_')[1:3]
			mkdir('/'.join([out_dir,name,ch])
			get_bcd_list(radecfile,data_dir,aors,ch,hdr)

	# now run map_to_sources to reverse the mapping such that each BCD has a 
	# list of its associated sources
	filepaths = glob.glob('bcd_dirs/*/*/bcd_list.json')
	pool.map(map_to_sources,filepaths)

	# now run get_bcd_phot to compute photometry on all the sources 
	source_list_paths = glob.glob('bcd_dirs/*/*/source_list.json')
	pool.map(get_bcd_phot,source_list_paths)

	# now run get_catalog to produce ch1/ch2,long/short catalogs
	phot_group_paths = glob.glob('bcd_dirs/*/*/phot_group.json')
	# NEW CODE HERE