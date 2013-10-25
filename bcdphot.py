#/usr/bin/env python
import os
import sys
import fnmatch
import multiprocessing
import glob
from get_bcd_list import get_bcd_list
from map_to_sources import map_bcd_to_sources
from get_bcd_phot import get_bcd_phot

def find_files(directory, pattern):
	for root, dirs, files in os.walk(directory):
		for basename in files:
			if fnmatch.fnmatch(basename, pattern):
				filename = os.path.join(root, basename)
				yield filename

if __name__ == "__main__":

	# instantiate the pool for parallel processing
	ncpus = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(processes=ncpus)
	print "using %i CPUs" % ncpus

	# this is the list of input sources (RA/Dec)
	infile = sys.argv[1]
	# this is the download directory containing all the data
	data_dir = sys.argv[2]

	# run get_bcd_list to associate individual BCDs to the input sources
	get_bcd_list(infile,data_dir)

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