#/usr/bin/env python

import os
import sys
import pyfits
import wcslib
import numpy as np
import simplejson as json
import multiprocessing
from scipy.spatial import cKDTree as KDT
from util import get_filepaths

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
	ra, dec, fitsfile = argslist
	hdr = pyfits.getheader(fitsfile)
	wcs = wcslib.WcsProjection(hdr)
	x, y = wcs.toPixel(ra, dec, rnd=False)
	x_max = hdr['NAXIS1']
	y_max = hdr['NAXIS2']
	if x >= 0 and x < x_max and y >= 0 and y < y_max:
		return (True,x,y)
	else:
		return (False,x,y)

def get_k_closest_bcd_idx(ra, dec, tree, k=10):
	coords = radec_to_coords(ra, dec)
	idx = tree.query(coords,k=k)[1].ravel()
	return idx

def radec_to_coords(ra, dec):
	x, y, z = spherical_to_cartesian(ra, dec)
	coords = np.empty((x.size, 3))
	coords[:, 0] = x
	coords[:, 1] = y
	coords[:, 2] = z
	return coords

def get_bcd_list(infile,data_dir,work_dir,aors,ch,hdr):

	# split the RA/Dec into two arrays
	radec = np.genfromtxt(infile)
	ra = radec[:,0]
	dec = radec[:,1]

	# populate the list of BCD files in the data dir
	# filenames = np.array([i for i in os.listdir(data_dir) if 'cbcd.fits' in i])
	# filepaths = np.array(['/'.join([data_dir,i]) for i in filenames])

	filepaths = np.array(get_filepaths('cbcd.fits',data_dir,aors,ch,hdr))
	filenames = np.array([i.split('/')[-1] for i in filepaths])

	# extract center pixel coordinates
	files_ra = []
	files_dec = []
	for fp in filepaths:
		hdr = pyfits.getheader(fp)
		center_ra = hdr['CRVAL1']
		center_dec = hdr['CRVAL2']
		files_ra.append(center_ra)
		files_dec.append(center_dec)

	# make arrays
	files_ra = np.array(files_ra, copy=False)
	files_dec = np.array(files_dec, copy=False)

	# make array of coordinates and populate the tree
	kdt = KDT(radec_to_coords(files_ra,files_dec))

	# spawn processes using multiprocessing to check for images containing,
	# the source, using the tree to find only the closest BCDs
	ncpus = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(processes=ncpus)
	# print "using %i CPUs" % ncpus

	max_num_images = 0
	sources = []
	for i in range(len(ra)):
		d = {'ra':ra[i], 'dec':dec[i]}
		print('finding files associated with source '+str(i+1)+\
			' at '+str(ra[i])+', '+str(dec[i]))
		# get the subset of BCDs to search
		idx = get_k_closest_bcd_idx(ra[i], dec[i], kdt, k=10)
		n_files = filepaths[idx].shape[0]
		filepaths_subset = filepaths[idx]
		filenames_subset = filenames[idx]
		argslist = zip([ra[i]]*n_files, [dec[i]]*n_files, filepaths_subset)
		# send jobs to the pool
		results = pool.map(source_in_image,argslist)
		# unzip the results and extract the boolean array and pixel coordinates
		results_unzipped = zip(*results)
		bool_arr = np.array(results_unzipped[0])
		x = results_unzipped[1]
		y = results_unzipped[2]
		pix_coord = np.array(zip(x,y))[bool_arr].tolist()
		# get the names of the files associated with the source
		good_bcds = filenames_subset[bool_arr].tolist()
		num_images = len(good_bcds)
		print('\t'+str(num_images)+' images')
		if num_images > max_num_images:
			max_num_images = num_images
		d['files'] = good_bcds
		d['pixels'] = pix_coord
		sources.append(d)

	outfile = 'bcd_list.json'
	outfilepath = '/'.join([work_dir,outfile])
	with open(outfilepath,'w') as w:
		json.dump(sources,w,indent=4*' ')
	print('created file: '+outfilepath)
	print('maximum number of images associated with a source: '+\
		str(max_num_images))


if __name__ == "__main__":

	# read the commandline arguments: list of RA/Dec and dir containing BCDs
	infile = sys.argv[1]
	data_dir = sys.argv[2]

	get_bcd_list(infile,data_dir)