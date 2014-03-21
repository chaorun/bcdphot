#/usr/bin/env python

import os
import sys
import pyfits
import pywcs
import itertools
import numpy as np
import simplejson as json
import multiprocessing
from scipy.spatial import cKDTree as KDT
from util import unzip, spherical_to_cartesian, radec_to_coords


def source_in_image(argslist):

	"""
	Reads the input FITS file and checks to make sure the RA/Dec pair
	corresponds to a physical pixel on the array.
	argslist = [RA, Dec, BCD]
	"""

	ra, dec, fitsfile = argslist
	hdr = pyfits.getheader(fitsfile)
	wcs = pywcs.WCS(hdr)
	x, y = map(lambda x: x[0], wcs.wcs_sky2pix(ra, dec, 1))
	x_max = hdr['NAXIS1']
	y_max = hdr['NAXIS2']
	if x >= 0 and x < x_max and y >= 0 and y < y_max:
		return True, x, y
	else:
		return False, x, y


def get_k_closest_bcd_idx(ra, dec, tree, k=10):

	"""
	Returns the indices of the k BCDs with central coordinates closest
	to the input RA/Dec coordinate pair.
	"""

	coords = radec_to_coords(ra, dec)
	idx = tree.query(coords, k=k)[1].ravel()
	return idx


def get_bcd_list(metadata):

	"""
	Metadata is a dict with keys:
		name, radecfile, data_dir, out_dir, work_dir, aors, channel,
		bcd_dict_path, max_cov
	"""

	radecfile = metadata['radecfile']
	work_dir = metadata['work_dir']
	aors = metadata['aors']
	max_cov = metadata['max_cov']

	# split the RA/Dec into two arrays
	radec = np.genfromtxt(radecfile)
	ra = radec[:,0]
	dec = radec[:,1]

	# read the region/ch/hdr specific bcd_dict in the work_dir for efficiency
	bcd_dict = json.load(open(metadata['bcd_dict_path']))
	filenames, filepaths = [np.array(i) for i in unzip(bcd_dict.items())]

	# extract center pixel coordinates
	files_ra = np.zeros(filepaths.size)
	files_dec = np.zeros(filepaths.size)
	for i, fp in enumerate(filepaths):
		hdr = pyfits.getheader(fp)
		files_ra[i] = hdr['CRVAL1']
		files_dec[i] = hdr['CRVAL2']

	# make array of coordinates and grow the tree
	kdt = KDT(radec_to_coords(files_ra, files_dec))

	# spawn processes using multiprocessing to check for images containing,
	# the source, using the tree to find only the closest BCDs to check
	ncpus = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(processes=ncpus)
	# print "using %i CPUs" % ncpus

	max_num_images = 0
	sources = []
	for i in range(len(ra)):

		# create internal source ID and associate with each RA/Dec pair
		d = {'id':i, 'ra':ra[i], 'dec':dec[i]}
		message = 'finding files associated with source {} at ({}, {})'
		print(message.format(i, ra[i], dec[i]))

		# get the subset of BCDs to search
		idx = get_k_closest_bcd_idx(ra[i], dec[i], kdt, k=max_cov)
		n_files = filepaths[idx].size
		filepaths_subset = filepaths[idx]
		filenames_subset = filenames[idx]
		argslist = zip([ra[i]]*n_files, [dec[i]]*n_files, filepaths_subset)

		# send jobs to the pool
		results = pool.map(source_in_image, argslist)

		# unzip the results and extract the boolean array and pixel coordinates
		results_unzipped = unzip(results)
		bool_arr = np.array(results_unzipped[0])
		
		# if none found, continue to next source
		if np.sum(bool_arr) == 0:
			continue

		x = results_unzipped[1]
		y = results_unzipped[2]
		pix_coord = np.array(zip(x, y))[bool_arr].tolist()

		# get the names of the files associated with the source
		good_bcds = filenames_subset[bool_arr].tolist()
		
		# compare the number of associated images to the previous maximum
		num_images = len(good_bcds)
		print('\t{} images'.format(num_images))
		if num_images > max_num_images:
			max_num_images = num_images

		# store results in source dict and append to source list
		d['files'] = good_bcds
		d['pixels'] = pix_coord
		sources.append(d)

	outfile = 'bcd_list.json'
	outfilepath = '/'.join([work_dir, outfile])
	with open(outfilepath, 'w') as w:
		json.dump(sources, w, indent=4*' ')

	print('created file: {}'.format(outfilepath))
	message = 'maximum number of images associated with a source: {}'
	print(message.format(max_num_images))


def map_bcd_sources(filepath):
	"""
	Reads JSON bcd list file, gets the set of BCD files,
	then associates each with a set of RA/Dec coordinates (and source ID).
	"""
	sources = json.load(open(filepath))

	list2d = [i['files'] for i in sources]
	merged = list(itertools.chain.from_iterable(list2d))
	bcd_list = list(set(merged))
	bcd_list.sort()

	d = {}
	for i in bcd_list:
		d[i] = []
		for s in sources:
			if i in s['files']:
				d[i].append((s['id'],s['ra'],s['dec']))

	outfilepath = filepath.replace('bcd_list.json','source_list.json')
	with open(outfilepath,'w') as w:
		json.dump(d,w,indent=4*' ')
	print('created file: '+outfilepath)
