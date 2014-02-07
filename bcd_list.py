#/usr/bin/env python

import os
import sys
import pyfits
import wcslib
import numpy as np
import simplejson as json
import multiprocessing
from scipy.spatial import cKDTree as KDT
from util import unzip

def spherical_to_cartesian(ra, dec):
	"""
	Inputs in degrees.  Outputs x,y,z
	"""
	rar = np.radians(ra)
	decr = np.radians(dec)
	x = np.cos(rar) * np.cos(decr)
	y = np.sin(rar) * np.cos(decr)
	z = np.sin(decr)
	return x, y, z

def source_in_image(argslist):
	"""
	Reads the input FITS file and checks to make sure the RA/Dec pair
	corresponds to a physical pixel on the array.
	argslist = [RA, Dec, BCD]
	"""
	ra, dec, fitsfile = argslist
	hdr = pyfits.getheader(fitsfile)
	wcs = wcslib.WcsProjection(hdr)
	x, y = wcs.toPixel(ra, dec, rnd=False)
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

def radec_to_coords(ra, dec):
	"""
	Helper function for constructing/querying k-d trees with coordinates
	in spherical geometry.
	Converts the input arrays from spherical coordinates to cartesion
	and populates a 3-dimensional array with the result.
	"""
	x, y, z = spherical_to_cartesian(ra, dec)
	coords = np.empty((x.size, 3))
	coords[:, 0] = x
	coords[:, 1] = y
	coords[:, 2] = z
	return coords

def get_bcd_list(metadata):
	"""
	Metadata is a dict with keys:
		name, radecfile, data_dir, out_dir, work_dir, aors, channel,
		bcd_dict_path, max_cov
	"""

	radecfile =	metadata['radecfile']
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