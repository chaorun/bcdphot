#/usr/bin/env python

import os
import sys
import glob
import pyfits
import numpy as np
import simplejson as json
import subprocess, shlex
from util import get_filepaths
from util import spherical_to_cartesian
from util import radec_to_coords
from util import great_circle_distance
from util import spherematch
from itertools import groupby


def get_k_closest_idx(ra, dec, tree, k=10):
	
	"""
	Converts RA/Dec to a cartesian coordinate array, then
	queries the input k-d tree to return nearest neighbors.
	Defaults to k=10 nearest neighbors.
	"""
	
	coords = radec_to_coords(ra, dec)
	idx = tree.query(coords,k=k)[1].ravel()
	return idx


def get_photometry_idl(source_list_path):
	
	"""
	Loop through the BCD files and get photometry on all associated sources.
	Returns arrays containing the output from bcd_phot.pro.
	"""

	work_dir = source_list_path.split('/source_list.json')[0]

	# read metadata for the working directory
	metadata = json.load(open(work_dir+'/metadata.json'))
	channel = metadata['channel']

	# set path to local IDL executable
	idl = metadata['idl_path']
	sources = json.load(open(source_list_path))
	bcd_dict = json.load(open(metadata['bcd_dict_path']))
	unc_dict = json.load(open(metadata['unc_dict_path']))
	msk_dict = json.load(open(metadata['msk_dict_path']))

	# initialize photometry output files with column names
	good_hdr = '# id ra dec ra_cen dec_cen x_cen y_cen flux_mjy unc_mjy '+\
		'flux_mjy_uncorrected'
	bad_hdr = '# id ra dec ra_cen dec_cen x_cen y_cen'
	with open(work_dir+'/good_list.txt','w') as g:
		g.write(good_hdr+'\n')
	with open(work_dir+'/bad_list.txt','w') as b:
		b.write(bad_hdr+'\n')
	if metadata['mask']:
		mask_hdr = '# id maskfile x y bitflag'
		with open(work_dir+'/masked_list.txt','w') as b:
			b.write(mask_hdr+'\n')

	# loop through the BCDs in the source list and get photometry
	for key in sources.keys():

		# keys are the BCD filenames
		bcd_path = bcd_dict[key]

		# get the UNC file path
		unc_key = key.replace('bcd.fits','bunc.fits')
		unc_path = unc_dict[unc_key]

		# also get the imask path
		msk_key = key.replace('cbcd.fits','bimsk.fits')
		msk_path = msk_dict[msk_key]

		# item of key is the list of ID/RA/Dec of sources in the image
		s = sources[key]

		# write to temp file so bcd_phot.pro can read it
		tmp_radec_path = work_dir+'/tmp_radec.txt'

		# np.savetxt(tmp_radec_path,s,fmt='%.9f')
		np.savetxt(tmp_radec_path,s,fmt=['%i']+['%.9f']*2)

		# spawn subprocess to get bcd_phot.pro output for the current image
		cmd = 'bcd_phot'+',"'+bcd_path+'","'+unc_path+'","'+msk_path+'","'+\
			tmp_radec_path+'",'+channel
		if metadata['mask']:
			cmd += ',/use_mask'
		returncode = subprocess.call([idl,'-quiet','-e',cmd], 
			stderr = subprocess.PIPE, stdout = subprocess.PIPE)

	# read the results of bcd_phot.pro
	print('created file: {}/good_list.txt'.format(work_dir))
	print('created file: {}/bad_list.txt'.format(work_dir))
	good_arr = np.loadtxt(work_dir+'/good_list.txt')
	bad_arr = np.loadtxt(work_dir+'/bad_list.txt')
	return good_arr, bad_arr


def get_phot_groups(good_arr):
	
	"""
	Use source ID numbers to collect photometry results of the same source
	into groups.
	"""

	# convert input array to list sorted by value of 'id' key of its elements
	good_lst_dct = sorted([{'id': int(i[0]), 'data': i[1:].tolist()} 
		for i in good_arr], key = lambda x: x['id'])

	# use a dict comprehension with groupby to get photometry groups
	phot_groups_dict = { key: [i['data'] for i in group] for 
		key, group in groupby(good_lst_dct, lambda x: x['id']) }
	return phot_groups_dict


def get_bcd_phot(source_list_path):
	
	"""
	Reads the output of map_bcd_sources() 'source_list.json' and 
	calls the functions get_photometry_idl() and get_phot_groups(),
	which in turn get photometry from IDL subprocesses and process
	the output.
	"""

	work_dir = source_list_path.split('/source_list.json')[0]

	# call get_photometry_idl() with output of map_bcd_sources()
	good_arr, bad_arr = get_photometry_idl(source_list_path)

	# collapse the output to their groupings, i.e. groups of 
	# measurements of the same star
	phot_groups_dict = get_phot_groups(good_arr)
	outfile = work_dir+'/phot_groups.json'
	with open(outfile,'w') as w:
		json.dump(phot_groups_dict,w,indent=4*' ')
	print('created file: '+outfile)

	# do the same for the failed measurements
	fail_groups_dict = get_phot_groups(bad_arr)
	outfile = work_dir+'/fail_groups.json'
	with open(outfile,'w') as w:
		json.dump(fail_groups_dict,w,indent=4*' ')
	print('created file: '+outfile)


def collapse_groups(phot_groups_dict):
	
	"""
	Computes the average RA, Dec, and flux, and the quadrature sum of the
	uncertainties for a group of photometric measurements of the same source,
	then returns a list of groups and their calculations.
	"""
	
	lst = []
	for key, value in phot_groups_dict.items():
		group = np.array(value)
		col_means = np.mean(group,0)

		# take the mean RA, Dec, and flux
		ra, dec = col_means[2:4]
		flux = col_means[6]
		flux_uncorrected = col_means[8]

		# sum the uncertainties in quadrature
		unc = np.sqrt(np.sum(group[:,7]**2))
		d = dict(id=key,ra=ra,dec=dec,flux=flux,unc=unc,group=value,
			flux_uncorrected=flux_uncorrected)
		lst.append(d)
	return lst


def write_mean_groups(phot_groups_path):
	
	"""
	Reads the output of get_bcd_phot() 'phot_groups.json' and collapses the
	groups of measurements to a single row per source
	"""
	
	phot_groups_dict = json.load(open(phot_groups_path))
	work_dir = phot_groups_path.split('/phot_groups.json')[0]
	phot_groups_mean = collapse_groups(phot_groups_dict)
	outfile = work_dir+'/phot_groups_mean.json'
	with open(outfile, 'w') as w:
		json.dump(phot_groups_mean, w, indent=' '*4)

def save_single_channel(phot_groups_mean_path):
	
	"""
	Creates a single channel/exposure catalog (no matching to other channel),
	and saves to disk.
	"""
	
	ch = json.load(open(phot_groups_mean_path))
	idnum = np.array( [ int(i['id']) for i in ch ] )
	ra = np.array( [ float(i['ra']) for i in ch ] )
	dec = np.array( [ float(i['dec']) for i in ch ] )
	flux = np.array( [ float(i['flux']) for i in ch ] )
	flux_uncorrected = np.array( [ float(i['flux_uncorrected']) for i in ch ] )
	unc = np.array( [ float(i['unc']) for i in ch ] )
	n_obs = np.array( [ len(i['group']) for i in ch ] )
	catalog = np.c_[idnum,ra,dec,flux,unc,flux_uncorrected,n_obs]
	work_dir = phot_groups_mean_path.split('/phot_groups_mean.json')[0]
	meta = json.load(open(work_dir+'/metadata.json'))
	if 'hdr' in meta.keys():
		out_name = '_'.join([meta['name'],meta['channel'],meta['hdr'],
			'catalog.txt'])
	else:
		out_name = '_'.join([meta['name'],meta['channel'],'catalog.txt'])
	out_path = '/'.join([work_dir,out_name])
	header = 'id ra dec flux unc flux_uncor n_obs'
	fmt = ['%i']+['%0.8f']*2+['%.4e']*3+['%i']
	idx = np.argsort(catalog[:,0])
	np.savetxt(out_path, catalog[idx], fmt = fmt, header = header)
	print("created file: "+out_path)


def save_catalog(catalog, out_path):
	
	"""
	Helper function for saving a matched catalog.
	"""

	header = 'ra dec ch1_flux ch1_unc '+\
		'ch2_flux ch2_unc n_obs1 n_obs2'
	np.savetxt(out_path, catalog, fmt = ['%.8f']*6+['%i']*2, header = header)
	print('created file: '+out_path)


def apply_array_location_correction(args_list):

	"""
	Matches the input ch1 and ch2 input (from write_mean_path()) and checks for
	sources with ch1>ch2, then applies the array location correction to these
	sources and saves the resulting matched ch1/ch2 catalog.
	"""

	ch1_path, ch2_path, out_path = args_list
	arrloc1 = pyfits.open('ch1_photcorr_ap_5.fits')[0].data
	arrloc2 = pyfits.open('ch2_photcorr_ap_5.fits')[0].data
	ch1 = np.array(json.load(open(ch1_path)))
	ch2 = np.array(json.load(open(ch2_path)))
	ra1 = np.array([i['ra'] for i in ch1])
	dec1 = np.array([i['dec'] for i in ch1])
	ra2 = np.array([i['ra'] for i in ch2])
	dec2 = np.array([i['dec'] for i in ch2])

	# match ch1/ch2 RA/Dec
	idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2, tolerance=1/3600.)
	ch1, ch2 = ch1[idx1], ch2[idx2]

	# get indices for the blue sources
	f1, f2 = [i['flux'] for i in ch1], [i['flux'] for i in ch2]

	# f1, f2 = [i['flux_uncorrected'] for i in ch1], [i['flux_uncorrected'] for i in ch2]
	blue = np.array(f1, copy=False) > np.array(f2, copy=False)

	# now loop through the matched sources and apply corrections
	catalog = []
	for i in range(ds.size):
		ra = np.mean( [ ch1[i]['ra'],ch2[i]['ra'] ] )
		dec = np.mean( [ ch1[i]['dec'],ch2[i]['dec'] ] )
		unc1, unc2 = ch1[i]['unc'], ch2[i]['unc']
		n_obs1, n_obs2 = len(ch1[i]['group']), len(ch2[i]['group'])
		is_blue = blue[i]
		if is_blue:
			group1 = np.array(ch1[i]['group'])
			coord1 = zip(group1[:,4], group1[:,5])
			group2 = np.array(ch2[i]['group'])
			coord2 = zip(group2[:,4], group2[:,5])
			for j in range(len(coord1)):
				x, y = [int(round(k)) for k in coord1[j]]
				group1[j,6:] *= arrloc1[x,y]
			for j in range(len(coord2)):
				x, y = [int(round(k)) for k in coord2[j]]
				group2[j,6:] *= arrloc2[x,y]
			flux1 = np.mean(group1[:,6])
			flux2 = np.mean(group2[:,6])
			unc1 = np.mean(group1[:,7])
			unc2 = np.mean(group2[:,7])
		else:
			flux1 = ch1[i]['flux']
			flux2 = ch2[i]['flux']
		row = [ra, dec, flux1, unc1, flux2, unc2, n_obs1, n_obs2]
		catalog.append(row)
	save_catalog(catalog, out_path)


def array_location_setup(filepaths):
	
	"""
	Takes a list of filepaths to the 'phot_groups_mean.json' files
	and organizes them by channel and exposure to create the list
	of arguments to be passed to apply_array_location_correction.
	"""
	
	ch1, ch2, out_paths = [], [], []
	for i in range(len(filepaths)):
		work_dir = filepaths[i].split('/phot_groups_mean.json')[0]
		metadata = json.load(open(work_dir+'/metadata.json'))
		channel = metadata['channel']
		if channel == '1':
			ch1.append(filepaths[i])
		elif channel == '2':
			ch2.append(filepaths[i])
	
	# make sure the ith elements of ch1 and ch2 are for the same region/exptime
	ch1.sort()
	ch2.sort()
	
	# construct file paths for matched ch1/ch2 catalogs
	for i in range(len(ch1)):
		work_dir1 = ch1[i].split('/phot_groups_mean.json')[0]
		meta1 = json.load(open(work_dir+'/metadata.json'))
		work_dir2 = ch2[i].split('/phot_groups_mean.json')[0]
		meta2 = json.load(open(work_dir+'/metadata.json'))
		out_dir = '/'.join( [ metadata['out_dir'], metadata['name'] ] )
		assert meta1['name'] == meta2['name'] 
		assert meta1['channel'] == meta2['channel']
		if 'hdr' in meta1.keys():
			assert meta1['hdr'] == meta2['hdr']
			out_name = '_'.join([meta1['name'],meta1['hdr'],
				'matched_catalog.txt'])
		else:
			out_name = '_'.join([meta1['name'],'matched_catalog.txt'])
		out_paths.append('/'.join([out_dir,out_name]))
	return zip(ch1, ch2, out_paths)