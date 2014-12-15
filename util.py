import os
import yaml
import fnmatch
import numpy as np
from functools import partial
from scipy.spatial import cKDTree as KDT


def find_files(directory, pattern):

	"""
	generator for filepaths in <directory> matching <pattern>
	"""

	for root, dirs, files in os.walk(directory):
		for basename in files:
			if fnmatch.fnmatch(basename, pattern):
				filename = os.path.join(root, basename)
				yield filename


def get_filepaths(suffix, data_dir, aors, ch, hdr=False):

	"""
	walks directory tree in <data_dir> and populates a list of filepaths
	matching the specifications given by <suffix>, <aors>, <ch>, <hdr>,
	where: 
		<suffix> : the filename ending to match, i.e. 'cbcd.fits'
		<aors> : a list of strings of AOR numbers
		<ch> : a string containing the channel number of IRAC 
		<hdr> : a string specifying the HDR exposure: 'long','short', or False
			if False, the data are assumed to be all the same exposure time
			(non-HDR mode)
	"""

	filepaths = []
	for aor in aors:
		for f in find_files(data_dir,'*I'+ch+'_'+aor+'*'+suffix):
			filename = f.split('/')[-1]
			expnum = filename.split('_')[3]
			if hdr == 'long':
				if int(expnum)%2 == 1:
					filepaths.append(f)
			elif hdr == 'short':
				if int(expnum)%2 == 0:
					filepaths.append(f)
			elif not hdr:
				filepaths.append(f)
			else:
				print('error: unexpected value for <hdr> argument')
	return filepaths


def ignore_oserror(f):

	"""
	decorator for handling OSError exceptions, i.e. when you try to
	make a directory that already exists.
	"""

	def wrapper(arg):
		try:
			f(arg)
		except OSError as e:
			pass
	return wrapper


@ignore_oserror
def mkdirs(dir):
	os.makedirs(dir)


def get_bcd_subset(bcd_dict, aors, ch, hdr=False):

	"""
	finds values (full paths) in bcd_dict with keys (file names) 
	matching aors, ch, and hdr
	"""

	bcds = []
	if not hdr:
		for key, value in bcd_dict.items():
			spl = key.split('_')
			if spl[1][1] == ch and spl[2] in aors:
				bcds.append(value)
		return bcds
	elif hdr == 'long':
		mod = 1
	elif hdr == 'short':
		mod = 0
	for key, value in bcd_dict.items():
		spl = key.split('_')
		if spl[1][1] == ch and spl[2] in aors and int(spl[3])%2 == mod:
			bcds.append(value)
	return bcds


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


def spherematch(ra1, dec1, ra2, dec2, tolerance=1/3600.):

	"""
	Uses a k-d tree to efficiently match two pairs of coordinates in spherical
	geometry, with a tolerance in degrees.
	"""

	args = ra1, dec1, ra2, dec2
	ra1, dec1, ra2, dec2 = map(partial(np.array, copy=False), args)
	coords1 = radec_to_coords(ra1, dec1)
	coords2 = radec_to_coords(ra2, dec2)
	kdt = KDT(coords2)
	idx2 = kdt.query(coords1)[1]
	ds = great_circle_distance(ra1, dec1, ra2[idx2], dec2[idx2])
	idx1 = np.arange(ra1.size)
	msk = ds < tolerance
	idx1 = idx1[msk]
	idx2 = idx2[msk]
	ds = ds[msk]
	return idx1, idx2, ds


def unzip(list):

	"""
	The inverse of zip()
	"""

	return zip(*list)


def discrete_hist(x):

	"""
	Histogram generator for discrete valued list input: yields value, count
	i.e. 
	for value, count in discrete_hist(x):
		print "{}: {}".format(value,count)
	"""

	s = set(x)
	arr = np.array(x)
	for i in s:
		yield i, np.sum(arr==i)


def to_mags(flux, zp):

	"""
	Converts flux to magnitudes using the input zero point,
	which must be in the same units as flux (i.e. Janskies)
	"""

	return -2.5 * np.log10(flux/zp)


def spz_jy_to_mags(jy, ch):

	"""
	Converts IRAC ch1 and ch2 flux in Janskies to magnitudes.
	"""

	if ch==1:
		zp = 280.9
	elif ch==2:
		zp = 179.7
	else:
		raise ValueError('ch must be either 1 or 2')
	return to_mags(jy,zp)


def ordinary_least_squares(y, X, intercept=False):

	"""
	Fits an ordinary least squares regression model. X is the
	design matrix and y is the response variable. If intercept
	equals True then an additional column of ones is added to the
	design matrix so that a non-zero intercept term is allowed.
	
	Returns an array of the parameters where the order corresponds
	to the columns in the design matrix (with the first being the
	intercept term if intercept=True), or a single float for the
	slope in the case of a 1-dimensional design matrix (vector)
	with intercept=False.
	"""

	if intercept:
		X = np.c_[np.ones(X.shape[0]), X]
	if len(X.shape) == 1:
		return np.dot( np.dot( 1. / ( np.dot(X.T, X) ), X.T), y )
	else:
		return np.dot( np.dot( np.linalg.inv( np.dot(X.T, X) ), X.T), y )


def setup_output_dirs(setup):

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

	if setup['params']['snr_dist_cull']:
		# get the minimum SNR value for an individual measurement to be kept
		min_snr = params['min_snr']
		# get the maximum distance from input RA/Dec to centroid position [arcsec]
		max_dist = params['max_dist']

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
	bcd_paths = list(find_files(data_dir, bcd_suffix))
	bcd_dict = {i.split('/')[-1]:i for i in bcd_paths}
	with open(out_dir+'/bcd_dict.json','w') as w:
		json.dump(bcd_dict, w, indent=' '*4)
	# do the same for UNC files
	unc_paths = list(find_files(data_dir, unc_suffix))
	unc_dict = {i.split('/')[-1]:i for i in unc_paths}
	with open(out_dir+'/unc_dict.json','w') as w:
		json.dump(unc_dict, w, indent=' '*4)
	# do the same for MSK files
	msk_paths = list(find_files(data_dir, '*_bimsk.fits'))
	msk_dict = {i.split('/')[-1]:i for i in msk_paths}
	with open(out_dir+'/msk_dict.json','w') as w:
		json.dump(msk_dict, w, indent=' '*4)

	# loop through the source lists (radecfiles) and create output directory
	# structure and metadata files used throughout the rest of the pipeline
	print('creating directory structure and metadata in output directory...')
	for region in regions:

		name = region['name']
		radecfiles = region['radec']
		aors = region['aors']

		for radec in radecfiles:
			f = radec['filepath']
			if is_hdr:
				ch = str(radec['channel'])
				hdr = radec['hdr']
				work_dir = '/'.join([out_dir, name, ch, hdr])
				mkdirs(work_dir)
			else:
				ch = str(radec['channel'])
				work_dir = '/'.join([out_dir, name, ch])
				mkdirs(work_dir)

			# make subset bcd_dict/unc_dict in each work_dir
			if is_hdr:
				bcd_paths_subset = get_bcd_subset(bcd_dict, aors, ch, hdr)
				unc_paths_subset = get_bcd_subset(unc_dict, aors, ch, hdr)
				msk_paths_subset = get_bcd_subset(msk_dict, aors, ch, hdr)
			else:
				bcd_paths_subset = get_bcd_subset(bcd_dict, aors, ch)
				unc_paths_subset = get_bcd_subset(unc_dict, aors, ch)
				msk_paths_subset = get_bcd_subset(msk_dict, aors, ch)

			bcd_dict_subset = {i.split('/')[-1]:i for i in bcd_paths_subset}
			unc_dict_subset = {i.split('/')[-1]:i for i in unc_paths_subset}
			msk_dict_subset = {i.split('/')[-1]:i for i in msk_paths_subset}

			bcd_dict_path = work_dir+'/bcd_dict.json'
			unc_dict_path = work_dir+'/unc_dict.json'
			msk_dict_path = work_dir+'/msk_dict.json'

			with open(bcd_dict_path,'w') as w:
				json.dump(bcd_dict_subset, w, indent=' '*4)
			print('created: '+bcd_dict_path)

			with open(unc_dict_path,'w') as w:
				json.dump(unc_dict_subset, w, indent=' '*4)
			print('created: '+unc_dict_path)

			with open(msk_dict_path,'w') as w:
				json.dump(msk_dict_subset, w, indent=' '*4)
			print('created: '+msk_dict_path)

			metadata = {
				'name': name,
				'data_dir':data_dir,
				'work_dir':work_dir,
				'out_dir': out_dir, 
				'radecfile': radecfile,
				'bcd_dict_path': bcd_dict_path,
				'unc_dict_path': unc_dict_path,
				'msk_dict_path': msk_dict_path,
				'aors':aors, 
				'channel':ch, 
				'cbcd': params['cbcd'],
				'mask': params['mask'],
				'idl_path': params['idl_path'], 
				'max_cov': params['max_cov'],
				'sigma_clip':params['sigma_clip'],
				'centroid':params['centroid']
				}
			if setup['params']['snr_dist_cull']:
				metadata['min_snr'] = min_snr
				metadata['max_dist'] = max_dist
			if is_hdr:
				if ch == '1':
					metadata['long_cutoff'] = params['long_cutoff_ch1']
					metadata['short_cutoff'] = params['short_cutoff_ch1']
				elif ch == '2':
					metadata['long_cutoff'] = params['long_cutoff_ch2']
					metadata['short_cutoff'] = params['short_cutoff_ch2']
				metadata['hdr'] = hdr
			metadata_path = work_dir+'/metadata.json'
			with open(metadata_path,'w') as w:
				json.dump(metadata, w, indent=' '*4)
			print('created: '+metadata_path)
