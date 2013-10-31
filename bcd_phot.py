#/usr/bin/env python

import os
import sys
import glob
import pyfits
import wcslib
import numpy as np
import simplejson as json
import subprocess, shlex
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

def great_circle_distance(ra1, dec1, ra2, dec2):
	"""
	Returns great circle distance.  Inputs in degrees.
	Uses vicenty distance formula - a bit slower than others, but
	numerically stable.
	"""
	from numpy import radians, degrees, sin, cos, arctan2, hypot
	# terminology from the Vicenty formula - lambda and phi and
	# "standpoint" and "forepoint"
	lambs = radians(ra1)
	phis = radians(dec1)
	lambf = radians(ra2)
	phif = radians(dec2)
	dlamb = lambf - lambs
	numera = cos(phif) * sin(dlamb)
	numerb = cos(phis) * sin(phif) - sin(phis) * cos(phif) * cos(dlamb)
	numer = hypot(numera, numerb)
	denom = sin(phis) * sin(phif) + cos(phis) * cos(phif) * cos(dlamb)
	return degrees(arctan2(numer, denom))

def get_k_closest_idx(ra, dec, tree, k=10):
	"""
	Converts RA/Dec to a cartesian coordinate array, then
	queries the input k-d tree to return nearest neighbors.
	Defaults to k=10 nearest neighbors.
	"""
	coords = radec_to_coords(ra, dec)
	idx = tree.query(coords,k=k)[1].ravel()
	return idx

def radec_to_coords(ra, dec):
	"""
	Converts the input RA/Dec from spherical degrees to cartesian,
	then returns an array containing the result.
	"""
	x, y, z = spherical_to_cartesian(ra, dec)
	# this is equivalent to, but faster than just doing np.array([x, y, z])
	coords = np.empty((x.size, 3))
	coords[:, 0] = x
	coords[:, 1] = y
	coords[:, 2] = z
	return coords	

def get_gross_list(source_list_path,metadata):
	"""
	Loop through the BCD files and get photometry on all associated sources.
	The result will be a 'gross' list of all the output from bcd_phot.pro.
	"""

	proj_dir,work_dir,bcd_paths,unc_paths,channel = 	metadata['proj_dir'],\
														metadata['work_dir'],\
														metadata['bcd_paths'],\
														metadata['unc_paths'],\
														metadata['channel']

	idl = '/usr/admin/local/itt/idl70/bin/idl'
	sources = json.load(open(source_list_path))
	# data_dir = source_list_path.split('source_list.json')[0]
	
	# create a dictionary with BCD filenames as the keys, full paths as values
	# bcd_paths = get_filepaths('cbcd.fits',proj_dir,'*','*')
	# bcd_dict = {i.split('/')[-1]:i for i in bcd_paths}

	# cbunc_paths = get_filepaths('cbunc.fits',proj_dir,'*','*')
	# cbunc_dict = {i.split('/')[-1]:i for i in cbunc_paths}

	# channel = source_list_path.split('_')[2][2]

	bcd_dict = {i.split('/')[-1]:i for i in bcd_paths}
	unc_dict = {i.split('/')[-1]:i for i in unc_paths}

	gross_lst = []
	for key in sources.keys():
		# keys are the BCD filenames
		bcd_path = bcd_dict[key]
		unc_key = key.replace('_cbcd.fits','_cbunc.fits')
		unc_path = unc_dict[unc_key]

		# now find the corresponding *cbunc.fits file to pass to bcd_phot.pro
		# for i in cbuncpaths:
		# 	if key_base in i:
		# 		unc_path = i
		# 		break
		# cbunc_path = [v for k,v in cbunc_dict.iteritems() \
		# 	if key_base in k][0]

		# item for key is the list of RA/Dec of sources in the image
		s = sources[key]
		# write to temp file so bcd_phot.pro can read it
		tmp_radec_path = work_dir+'/tmp_radec.txt'
		np.savetxt(tmp_radec_path,s,fmt='%.9f')
		# spawn subprocess to get bcd_phot.pro output for the current image
		cmd = 'echo bcd_phot,"'+bcd_path+'","'+unc_path+'","'+tmp_radec_path+\
			'",'+channel
		p1 = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
		p2 = subprocess.Popen([idl], stdin=p1.stdout, stdout=subprocess.PIPE)
		# get the bcd_phot.pro output as a single string
		result_str = p2.stdout.read()
		# split the string on newlines
		result_split = result_str.strip().split('\n')
		# split each line on whitespace
		result_lst = [i.split() for i in result_split if 'NaN' not in i]
		# append the result list to the gross list
		gross_lst = gross_lst+result_lst
	return gross_lst

def get_phot_groups(gross_arr):
	""" 
	Uses a k-d tree to get groupings of measurements from gross array.
	A tolerance of 1 arcsec is the criterion for group membership used
	to produce the 'net' list, which contains arrays of the groups.
	"""
	# construct the k-d tree from the sources' cartesian coordinates
	ra, dec = gross_arr[:,0],gross_arr[:,1]
	coords = radec_to_coords(ra, dec)
	kdt = KDT(coords)
	# set tolerance to 1 arcsec (in degrees)
	tol = 1/3600.
	# list to store the groups of bcd_phot.pro output
	phot_groups = []
	# list to store the hashes of each sorted group array
	hashes = []
	for i in range(gross_arr.shape[0]):
		# get the ith pair of RA/Dec to search the k-d tree with
		rai, deci = gross_arr[i,0],gross_arr[i,1]
		# search the tree
		idx = get_k_closest_idx(rai,deci,kdt,k=10)
		# compute distances of query result RA/Dec to search RA/Dec
		ds = great_circle_distance(rai, deci, ra[idx], dec[idx])
		# mask out the query results more distant than the tolerance
		msk = ds < tol
		idx = idx[msk]
		# sort the index array so order is preserved in each query result
		#	(this allows the hashing trick to work)
		idx.sort()
		# compute hash of string reprentation of result array (fastest)
		h = hash(gross_arr[idx,:].tostring())
		# only keep the results (measurement groupings) we haven't seen already
		if h not in hashes:
			hashes.append(h)
			phot_groups.append(gross_arr[idx,:])
	phot_groups_lists = [i.tolist() for i in phot_groups]
	phot_groups_dict = dict(zip(hashes,phot_groups_lists))
	return phot_groups_dict

def get_bcd_phot(source_list_path):
	"""
	Reads the output of map_bcd_sources() 'source_list.json' and 
	calls the functions get_gross_list() and get_phot_groups(),
	which in turn get photometry from IDL subprocesses and process
	the output.
	"""
	work_dir = source_list_path.split('/source_list.json')[0]
	# read metadata for the working directory
	metadata = json.load(open(work_dir+'/metadata.json'))
	# call get_gross_list() with output of map_bcd_sources() and metadata
	gross_lst = get_gross_list(source_list_path,metadata)
	# save the gross output (from IDL) to text
	gross_arr = np.array(gross_lst).astype(np.float)
	outfile = work_dir+'/gross_arr.txt'
	np.savetxt(outfile,gross_arr,fmt='%.18e')
	print('created file: '+outfile)
	# collapse the gross output to their groupings, i.e. groups of 
	# measurements of the same star
	phot_groups_dict = get_phot_groups(gross_arr)
	outfile = work_dir+'/phot_groups.json'
	with open(outfile,'w') as w:
		json.dump(phot_groups_dict,w,indent=4*' ')
	print('created file: '+outfile)

def collapse_groups(phot_groups_dict):
	"""
	Computes the average RA, Dec, and flux, and the quadrature sum of the
	uncertainties for a group of photometric measurements of the same source.
	"""
	lst = []
	for value in phot_groups_dict.values():
		group = np.array(value)
		col_means = np.mean(group,0)
		ra, dec = col_means[:2]
		flux = col_means[4]
		unc = np.sqrt(np.sum(group[:,5]**2))
		row = [ra,dec,flux,unc]
		lst.append(row)
	return lst

def save_catalog(phot_groups_path):
	"""
	Reads the output of get_bcd_phot() 'phot_groups.json' and collapses the
	groups of measurements to a single row per source
	"""
	phot_groups_dict = json.load(open(phot_groups_path))
	work_dir = phot_groups_path.split('/phot_groups.json')[0]
	phot_arr = np.array(collapse_groups(phot_groups_dict))
	outfile = work_dir+'/phot_collapsed.txt'
	np.savetxt(outfile,phot_arr,fmt='%.18e')
	# return phot_arr

def apply_array_location_correction(ch1_cat_path,ch2_cat_path):
	ch1_cat = np.genfromtxt(ch1_cat_path)
	ch2_cat = np.genfromtxt(ch2_cat_path)
	ch1_flux = ch1_cat[:,2]
	ch2_flux = ch2_cat[:,2]
	# need to match ch1 and ch2 sources by RA/Dec