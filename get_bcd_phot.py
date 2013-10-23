import os
import sys
import glob
import pyfits
sys.path.append('/home/jlivings/python')
# sys.path.append('/Users/jlivings/python_lib')
import wcslib
import numpy as np
import simplejson as json
import subprocess, shlex
from scipy.spatial import cKDTree as KDT

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

def get_gross_list(source_list_path,cbuncpaths):
	"""
	Loop through the BCD files and get photometry on all associated sources.
	The result will be a 'gross' list of all the output from bcd_phot.pro.
	"""
	idl = '/usr/admin/local/itt/idl70/bin/idl'
	sources = json.load(open(source_list_path))
	data_dir = source_list_path.split('source_list.json')[0]
	channel = source_list_path.split('_')[2][2]
	gross_lst = []
	for key in sources.keys():
		# keys are the BCD filenames
		bcd_path = data_dir+key
		key_base = key.split('_cbcd.fits')[0]
		# now find the corresponding *cbunc.fits file to pass to bcd_phot.pro
		for i in cbuncpaths:
			if key_base in i:
				unc_path = i
				break
		# item for key is the list of RA/Dec of sources in the image
		s = sources[key]
		# write to temp file so bcd_phot.pro can read it
		tmp_radec_path = data_dir+'tmp_radec.txt'
		np.savetxt(tmp_radec_path,s,fmt='%.9f')
		# spawn subprocess to get bcd_phot.pro output for the current image
		cmd = 'echo bcd_phot,"'+bcd_path+'","'+unc_path+'","'+tmp_radec_path+
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

def get_phot_groups(gross_arr,data_dir):
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
	d = dict(zip(hashes,phot_groups_lists))
	json.dump(d,open(data_dir+'phot_group.json','w'),indent=4*' ')
	print('created file: '+data_dir+'phot_group.json')
	return phot_groups


if __name__ == "__main__":
	source_list_path = sys.argv[1]
	# source_list = 'bcd_dirs/d765_ch1/short/source_list.json'
	data_dir = source_list_path.split('source_list.json')[0]
	idl = '/usr/admin/local/itt/idl70/bin/idl'

	cbuncpaths = glob.glob('unzipdirs/*/*/*/*cbunc.fits')

	gross_lst = get_gross_list(source_list_path,cbuncpaths)
	# make an array from the gross list and save to data_dir
	gross_arr = np.array(gross_lst).astype(np.float)
	np.savetxt(data_dir+'gross_arr.txt',gross_arr)
	print('created file: '+data_dir+'gross_arr.txt')
	phot_groups = get_phot_groups(gross_arr,data_dir)
