import os
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


def get_filepaths(suffix,data_dir,aors,ch,hdr=False):

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


def get_bcd_subset(bcd_dict,aors,ch,hdr=False):

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


