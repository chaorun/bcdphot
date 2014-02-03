import os
import fnmatch
import numpy as np

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

# unnecessary, just for fun
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
