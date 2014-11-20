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
from util import ordinary_least_squares
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
		if metadata['centroid']:
			cmd += ',/centroid'
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


def cull_bad_measurements(phot_groups_filepath):

	"""
	Eliminates bad measurements from the groups of measurements
	for each source. Uses a SNR cutoff and a proximity cutoff, 
	such that only sources with SNR above min_snr and with centroids
	less than max_dist from the input RA/Dec position will survive.
	max_dist is in units of arcsec.
	"""

	work_dir = '/'.join(phot_groups_filepath.split('/')[:-1])
	meta = json.load(open(work_dir+'/metadata.json'))

	min_snr = meta['min_snr']
	max_dist = meta['max_dist']

	# read in the photometry JSON files
	ch = json.load(open(phot_groups_filepath))

	# loop through the sources and look for measurements to cull
	rejected = {}
	badkeys = []
	for key in ch:
		good, bad = [], []
		for obs in ch[key]:
			ra_inp, dec_inp = obs[:2]
			ra_cnt, dec_cnt = obs[2:4]
			snr = obs[6]/obs[7]
			d = great_circle_distance(ra_inp, dec_inp, ra_cnt, dec_cnt) * 3600
			if snr > min_snr and d < max_dist:
				good.append(obs)
			else:
				bad.append(obs)
		if len(good) > 0:
			ch[key] = good
		else:
			badkeys.append(key)
		if len(bad) > 0:
			rejected[key] = bad
	# eliminate source altogether if all its measurements get culled
	for key in badkeys:
		ch.pop(key)

	# write to disk
	out_path = work_dir+'/phot_groups_culled.json'
	with open(out_path,'w') as w:
		json.dump(ch, w, indent=4*' ')
	print('created file: '+out_path)
	out_path = work_dir+'/phot_groups_rejected.json'
	with open(out_path,'w') as w:
		json.dump(rejected, w, indent=4*' ')
	print('created file: '+out_path)


def apply_array_location_correction(phot_groups_filepath):

	"""
	Reads the phot_groups.json files and applies the array 
	location correction to all measurements of ALL sources.
	Writes the result to disk in work_dir as
	'phot_groups_arrayloc.json'
	"""

	work_dir = '/'.join(phot_groups_filepath.split('/')[:-1])
	meta = json.load(open(work_dir+'/metadata.json'))

	# read in the array location correction values
	if meta['channel'] is '1':
		arrloc = pyfits.open('ch1_photcorr_ap_5.fits')[0].data
	elif meta['channel'] is '2':
		arrloc = pyfits.open('ch2_photcorr_ap_5.fits')[0].data

	# read in the photometry JSON files
	ch = json.load(open(phot_groups_filepath))

	# apply correction
	for key in ch:
		for obs in ch[key]:
			x, y = obs[4:6]
			obs[6:] = [i * arrloc[x,y] for i in obs[6:]]

	# write to disk
	out_path = work_dir+'/phot_groups_arrayloc.json'
	with open(out_path,'w') as w:
		json.dump(ch, w, indent=4*' ')
	print('created file: '+out_path)


def uncorrect_red_sources(phot_groups_filepath_tuple):

	"""	
	identify the 'red' sources by comparing ch1 to ch2 energies, then find their
	entries in the single-exposure catalogs and un-correct (divide) the
	array location dependent correction. red = ch1_flux < 3.6/4.5 * ch2_flux
	"""	

	ch1_file, ch2_file = phot_groups_filepath_tuple

	arrloc1 = pyfits.open('ch1_photcorr_ap_5.fits')[0].data
	arrloc2 = pyfits.open('ch2_photcorr_ap_5.fits')[0].data

	ch1 = json.load(open(ch1_file))
	ch2 = json.load(open(ch2_file))

	ra1, dec1 = zip(*[i[0][:2] for i in ch1.values()])
	ra2, dec2 = zip(*[i[0][:2] for i in ch2.values()])

	if len(ra1) < len(ra2):
		idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2, tolerance=2/3600.)
	else:
		idx2, idx1, ds = spherematch(ra2, dec2, ra1, dec1, tolerance=2/3600.)

	for k1, k2 in zip(np.array(ch1.keys())[idx1], np.array(ch2.keys())[idx2]):
		f1 = np.array(ch1[k1]).mean(0)[6]
		f2 = np.array(ch2[k2]).mean(0)[6]
		red = f1 < (3.6/4.5) * f2
		if red:
			for obs in ch1[k1]:
				x, y = obs[4:6]
				obs[6:] = [i / arrloc1[x,y] for i in obs[6:]]
			for obs in ch2[k2]:
				x, y = obs[4:6]
				obs[6:] = [i / arrloc2[x,y] for i in obs[6:]]

	out_path = ch1_file.replace('_arrayloc.json', '_arrayloc_cor.json')
	with open(out_path,'w') as w:
		json.dump(ch1, w, indent=4*' ')
	print('created file: '+out_path)
	out_path = ch2_file.replace('_arrayloc.json', '_arrayloc_cor.json')
	with open(out_path,'w') as w:
		json.dump(ch2, w, indent=4*' ')
	print('created file: '+out_path)


def calculate_full_uncertainties(phot_groups_filepath):

	"""
	Reads phot_groups.json or phot_groups_arrayloc.json files
	and calculates the full uncertainties including both systematic
	and photometric noise sources. The calculation is as follows:

		sigma_phot = sqrt(sum(deltaF_i**2)) / n_obs
		sigma_sys = median MAD of 99th percentile brightest flux w/ n_obs > 2
		sigma_tot = sqrt(sigma_phot**2 + sigma_sys**2)

	where MAD = median absolute deviation = median(abs(x-median(x)))

	Outputs a catalog with the mean flux and full uncertainty values
	"""

	work_dir = '/'.join(phot_groups_filepath.split('/')[:-1])
	meta = json.load(open(work_dir+'/metadata.json'))

	# read in the photometry JSON files
	ch = json.load(open(phot_groups_filepath))

	# calculate the systematic uncertainties
	flux, n_obs, mad = [], [], []
	for key in ch:
		flux_i = [obs[6] for obs in ch[key]]
		flux.append( np.mean( flux_i ) )
		n_obs.append( len(ch[key]) )
		mad.append( np.median(np.abs(flux_i-np.median(flux_i))) )
	flux, n_obs, mad = map(np.array, (flux, n_obs, mad))
	flux_n2 = flux[n_obs > 2]
	idx = np.argsort(flux_n2)
	min_flux = flux_n2[idx][-100]
	brightest = (flux >= min_flux) & (n_obs > 2)
	# brightest = (flux > np.percentile(flux, 99)) & (n_obs > 1)
	sigma_sys = np.median(mad[brightest] / flux[brightest])

	# print the systematic uncertainty to stdout
	if 'hdr' in meta.keys():
		msg = "region: {}, channel: {}, exposure: {}\n"+\
		"systematic uncertainty: {}\n"+\
		"calculated from {} datapoints\n"+\
		"average number of measurements: {}"
		msg = msg.format(meta['name'], meta['channel'], meta['hdr'],
			sigma_sys, brightest.sum(), n_obs[brightest].mean())
		with open(work_dir+'/systematic_uncertainty.txt','w') as w:
			w.write(msg)
	else:
		msg = "region: {}, channel: {}\n"+\
		"systematic uncertainty: {}\n"+\
		"calculated from {} datapoints\n"+\
		"average number of measurements: {}"
		msg = msg.format(meta['name'], meta['channel'], 
			sigma_sys, brightest.sum(), n_obs[brightest].mean())
		with open(work_dir+'/systematic_uncertainty.txt','w') as w:
			w.write(msg)

	# calculate the photometric uncertainties
	sigma_phot = []
	for key in ch:
		phot_unc = np.array( [obs[7]/obs[6] for obs in ch[key]] )
		sigma_phot.append( np.sqrt(np.sum(phot_unc**2)) / phot_unc.size )
	sigma_phot = np.array(sigma_phot)

	# calculate full uncertainties
	sigma_tot = np.sqrt(sigma_phot**2 + sigma_sys**2)

	# write to disk
	# columns: id, ra, dec, flux, unc, n_obs
	ids, ra, dec = [], [], []
	for key in ch:
		ids.append(int(key))
		ra.append(np.mean([obs[2] for obs in ch[key]]))
		dec.append(np.mean([obs[3] for obs in ch[key]]))
	data = np.c_[ids, ra, dec, flux, flux*sigma_tot, n_obs]
	header = 'id ra dec flux unc n_obs'
	fmt = ['%i']+['%0.8f']*2+['%.4e']*2+['%i']	
	idx = np.argsort(ids)
	if 'hdr' in meta.keys():
		out_name = '_'.join([meta['name'], meta['channel'], meta['hdr'],
			'catalog.txt'])
	else:
		out_name = '_'.join([meta['name'], meta['channel'], 'catalog.txt'])
	out_path = '/'.join([work_dir, out_name])
	np.savetxt(out_path, data[idx], fmt = fmt, header = header)
	print('created file: '+out_path)


def combine_hdr_catalogs(catalog_filepaths_tuple):

	"""
	Takes a tuple containing the filepaths to the short and long exposure
	single-channel catalogs for a given region and channel. The result is
	a single catalog containing the union of all sources in both short and
	long exposure catalogs, with the short exposure measurements being used
	for the brighter sources, and the long exposure measurements used for 
	the fainter sources. The cutoff between the two is determined by the 
	parameter 'hdr_cutoff' in the metadata file and should be set to the 
	saturation limit for the long exposure data (there are actually 2 
	parameters, one for each channel: hdr_cutoff_ch1, hdr_cutoff_ch2).
	"""

	# read in the data
	long_file, short_file = catalog_filepaths_tuple
	work_dir = '/'.join(short_file.split('/')[:-1])
	meta = json.load(open(work_dir+'/metadata.json'))
	header = 'id ra dec flux unc n_obs'
	names = header.split()
	long_cat = np.recfromtxt(long_file, names=names)
	short_cat = np.recfromtxt(short_file, names=names)

	# fit a line to short ~ long
	idx_s = short_cat.flux < meta['short_cutoff']
	idx_l = long_cat.flux < meta['long_cutoff']
	short_flux = short_cat.flux[idx_s]
	long_flux = long_cat.flux[idx_l]
	short_ra, short_dec = short_cat.ra[idx_s], short_cat.dec[idx_s]
	long_ra, long_dec = long_cat.ra[idx_l], long_cat.dec[idx_l]
	idx1, idx2, ds = spherematch(short_ra, short_dec, long_ra, 
		long_dec, tolerance=1/3600.)
	y = short_flux[idx1]
	X = long_flux[idx2]
	slope = ordinary_least_squares(y, X)

	# divide short flux/unc by the slope so that it agrees with the long flux
	print('region {} correction value: {}'.format(meta['name'], slope))
	short_cat.flux /= slope
	short_cat.unc /= slope

	# get everything brighter than the cutoff in short and combine with long
	idx_faint = long_cat.flux < meta['long_cutoff']
	idx_bright = short_cat.flux > meta['long_cutoff']

	# before concatenation of long and short subsets, check for any duplicates
	# (if they exist they should tend to have flux very close to the cutoff)
	ls, ss = long_cat[idx_faint], short_cat[idx_bright]
	idx_s, idx_l, ds = spherematch(ss.ra, ss.dec, ls.ra, ls.dec, 
		tolerance=1/3600.)
	dup_ids = []
	for idx in idx_l:
		if (ls.flux[idx] > 0.9 * meta['long_cutoff']) & \
			(ls.flux[idx] < meta['long_cutoff']):
			dup_ids.append(ls.id[idx])
	
	# now use the ids of the duplicates to delete them from the long dataset
	for idx in dup_ids:
		ls = ls[ls.id != idx]

	data = np.concatenate([ls, ss])

	# eliminate sources with negative flux
	good = data['flux'] >= 0
	data = data[good]

	# apply global sigma clip using the value from setup.yaml
	snr = data['flux'] / data['unc']
	good = snr >= meta['sigma_clip']
	data = data[good]

	# write to disk
	header = 'id ra dec flux unc n_obs'
	data = data[header.split()]
	idx = np.argsort(data['ra'])
	data = data[idx]
	data['id'] = np.arange(1, data.shape[0]+1)
	fmt = ['%i']+['%0.8f']*2+['%.4e']*2+['%i']
	out_name = '_'.join([meta['name'], meta['channel'], 
		'combined_hdr_catalog.txt'])
	out_path = '/'.join(['/'.join(work_dir.split('/')[:-1]), out_name])
	np.savetxt(out_path, data, fmt = fmt, header = header)
	print('created file: '+out_path)


def sigma_clip_non_hdr(filepath):

	"""
	Eliminates sources with SNR less than the 'sigma_clip' parameter from setup
	file. Also checks for negative flux sources which may remain after this
	step (can happen when uncertainty values are also negative).
	"""

	work_dir = '/'.join(filepath.split('/')[:-1])
	meta = json.load(open(work_dir+'/metadata.json'))

	names = open(filepath).readline().split()[1:]
	data = np.recfromtxt(filepath, names=names)

	# get rid of low SNR sources
	snr = data['flux'] / data['unc']
	good = snr >= meta['sigma_clip']
	data = data[good]

	# get rid of any remaining negative flux sources
	good = data['flux'] > 0
	data = data[good]

	# get rid of id column
	data = data[['ra', 'dec', 'flux', 'unc', 'n_obs']]

	# write to disk
	fmt = ['%0.8f']*2+['%.4e']*2+['%i']
	out_path = filepath.replace('.txt', '_sigclip.txt')
	header = ' '.join(names[1:])
	np.savetxt(out_path, data, fmt = fmt, header = header, comments='')
	print('created file: '+out_path)
