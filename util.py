import os
import fnmatch

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
