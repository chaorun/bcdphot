import os
import sys
import csv
import json
import numpy as np
import urllib2
import subprocess, shlex
from compiler.ast import flatten
from prettyplotlib import plt
from util import find_files
from util import spherical_to_cartesian
from util import radec_to_coords
from util import great_circle_distance
from scipy.spatial import cKDTree as KDT
from collections import OrderedDict


def format_radec(ra, dec):
	return "{}+{}".format(ra, dec)

def query_2mass_xsc_polygon(ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4):
	base_url = "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?outfmt=1"	
	query_url = base_url + "&spatial=polygon"
	query_url += "&polygon="+format_radec(ra1, dec1)+","
	query_url += format_radec(ra2, dec2)+","
	query_url += format_radec(ra3, dec3)+","
	query_url += format_radec(ra4, dec4)
	query_url += "&catalog=fp_xsc"
	return query_url

def get_region_corners(catalog):
	cat = np.recfromtxt(catalog, names=open(catalog).readline().split()[1:])
	c1 = (cat.ra.min(), cat.dec[cat.ra==cat.ra.min()][0])
	c2 = (cat.ra[cat.dec==cat.dec.max()][0], cat.dec.max())
	c3 = (cat.ra.max(), cat.dec[cat.ra==cat.ra.max()][0])
	c4 = (cat.ra[cat.dec==cat.dec.min()][0], cat.dec.min())
	return reduce(lambda x,y: x+y, [c1, c2, c3, c4])

def get_compliment(idx, size):
	return np.where(np.in1d(np.arange(size), idx, invert=True))[0]

def spherematch2(ra1, dec1, ra2, dec2, kdt, tolerance=1/3600., k=100):
	"""
	Uses a k-d tree to efficiently match two pairs of coordinates in spherical
	geometry, with a tolerance in degrees.
	"""
	coords1 = radec_to_coords(ra1, dec1)
	idx = kdt.query(coords1, k=k)[1].flatten()
	ds = great_circle_distance(ra1, dec1, ra2[idx], dec2[idx])
	msk = ds < tolerance
	idx = idx[msk]
	ds = ds[msk]
	return idx, ds

def check_n_in_aper(radius_factor=1, k=100):

	for catfile in find_files(bcdphot_out_path, "*_combined_hdr_catalog.txt"):

		print
		print catfile
		names = open(catfile).readline().split()[1:]
		cat = np.recfromtxt(catfile, names=names)

		xscfile = catfile.replace('combined_hdr_catalog.txt','2mass_xsc.tbl')
		print xscfile
		names = open(xscfile).read().split('\n')[76].split('|')[1:-1]
		xsc = np.recfromtxt(xscfile, skip_header=80, names=names)

		n_in_aper = []
		coords = radec_to_coords(cat.ra, cat.dec)
		kdt = KDT(coords)
		for i in range(xsc.size):
			r_deg = xsc.r_ext[i]/3600.

			idx, ds = spherematch2(xsc.ra[i], xsc.dec[i], cat.ra, cat.dec,
				kdt, tolerance=radius_factor*r_deg, k=k)
			n_in_aper.append(ds.size)
		for i in [(i,n_in_aper.count(i)) for i in set(n_in_aper)]:
			print i

# check_n_in_aper(k=500)
# --> the maximum k value required to handle all catalogs is 446


def run_xsc_phot(bcdphot_out_path, mosaic_path):
	replaced = {}
	for cat in find_files(bcdphot_out_path, "*_combined_hdr_catalog.txt"):

		print("\n======================================================")
		print("\nadjusting photometry in: {}".format(cat.split('/')[-1]))
		print("------------------------------------------------------")
		outpath = cat.replace('combined_hdr_catalog.txt','2mass_xsc.tbl')

		# retrieve 2mass data if file doesn't already exist (from previous run)
		if not os.path.isfile(outpath):
			# get url and retrieve data
			url = query_2mass_xsc_polygon(*get_region_corners(cat))
			print("\ndownloading 2MASS photometry from: {}".format(url))
			text = urllib2.urlopen(url).read()
			# write to disk
			with open(outpath, 'w') as f:
				f.write(text)
			print("\ncreated file: {}".format(outpath))

		# read back in as recarray	
		print("\nreading: {}".format(outpath))
		names = open(outpath).read().split('\n')[76].split('|')[1:-1]
		da = np.recfromtxt(outpath, skip_header=80, names=names)

		# write input file for xsc_phot.pro
		infile_outpath = '/'.join(cat.split('/')[:-1])+'/xsc.txt'
		with open(infile_outpath,'w') as w:
			for i in range(da.shape[0]):
				w.write("{} {} {} {}\n".format(da.designation[i], da.ra[i], da.dec[i], da.r_ext[i]))
		print("\ncreated input file for xsc_phot.pro: {}".format(infile_outpath))

		# locate the FITS mosaic file for xsc_phot.pro to do photometry on
		reg, ch = cat.split('/')[-1].split('_')[:2]
		mosaicfile = filter(lambda x: 'dirbe{}/ch{}/long/full/Combine'\
			.format(reg,ch) in x, find_files(mosaic_path, '*mosaic.fits'))[0]
		print("\nfound mosaic file: {}".format(mosaicfile))

		# spawn IDL subprocess running xsc_phot.pro and catch stdout in file
		outpath = infile_outpath.replace('xsc.txt', 'xsc_phot_out.txt')
		if not os.path.isfile(outpath):
			outfile = open(outpath,'w')
			print("\nspawning xsc_phot.pro IDL subprocess")
			cmd = "xsc_phot,'"+mosaicfile+"','"+infile_outpath+"','long'"
			rc = subprocess.call(['/usr/local/itt/idl71/bin/idl','-quiet','-e',cmd], 
				stderr = subprocess.PIPE, stdout = outfile)
			outfile.close()

		# read in output to recarray
		print("\nreading: {}".format(outpath))
		phot = np.recfromtxt(outpath, names=['id','flux','unc','sky','skyunc'])

		# make sure rows are aligned
		assert (da.designation == phot.id).all()

		# ignore xsc sources we got a NaN flux for
		bad = np.isnan(phot.flux)
		print("\naper.pro returned NaN flux for {} sources".format(bad.sum()))
		if bad.sum() > 0:
			for i in phot[bad].id:
				print(i)
			outpath = cat.replace('combined_hdr_catalog.txt','xsc_nan_phot.csv')
			with open(outpath,'w') as f:
				w = csv.writer(f)
				w.writerow(da.dtype.names)
				w.writerows(da[bad].tolist())
			print('\ncreated file: {}'.format(outpath))
		phot = phot[~bad]
		da = da[~bad]

		# read in pipeline catalog
		print("\nreading: {}".format(cat))
		names = open(cat).readline().split()[1:]
		c = np.recfromtxt(cat, names=names)

		# loop through xsc sources and find matches in pipeline catalog
		print("\nfinding records associated with XSC sources in pipeline catalog")
		c_flux_total = []
		n_in_aper = []
		c_idx = []
		coords = radec_to_coords(c.ra, c.dec)
		kdt = KDT(coords)
		for i in range(phot.size):
			radius = da.r_ext[i]/3600.
			# idx1, idx2, ds = spherematch(da.ra[i], da.dec[i], 
			# 	c.ra, c.dec, tolerance=radius)
			idx, ds = spherematch2(da.ra[i], da.dec[i], c.ra, c.dec,
				kdt, tolerance=radius, k=500)
			# c_flux_total.append(c.flux[idx2].sum())
			# n_in_aper.append(c.flux[idx2].size)
			# c_idx.append(idx2.tolist())
			c_flux_total.append(c.flux[idx].sum())
			n_in_aper.append(ds.size)
			c_idx.append(idx.tolist())
		print("\nhistogram of source counts in r_ext aperture")
		for i in [(i,n_in_aper.count(i)) for i in set(n_in_aper)]:
			print i

		# create new version of catalog file with xsc-associated entries replaced
		c_idx = np.array(flatten(c_idx))
		print("\nremoving {}, adding {}".format(c_idx.size, phot.size))
		replaced[cat] = {'old':c_idx.size, 'new':phot.size}
		replaced[cat]['hist'] = [(i,n_in_aper.count(i)) for i in set(n_in_aper)]
		c = np.delete(c, c_idx)
		newrows = np.rec.array([(-i, da.ra[i], da.dec[i], 
			phot.flux[i], phot.unc[i], 1) for i in \
			range(phot.size)], dtype=c.dtype)
		newcat = np.hstack((c, newrows))

		# write new version of catalog to disk
		fmt = ['%i']+['%0.8f']*2+['%.4e']*2+['%i']
		outpath = cat.replace('catalog.txt', 'catalog_xsc_cor.txt')
		np.savetxt(outpath, newcat, fmt = fmt, header = ' '.join(names))
		print('\ncreated file: {}'.format(outpath))

		# make plot of total old vs. new flux
		plt.scatter(c_flux_total, phot.flux)
		ylim = plt.gca().get_ylim()
		plt.xlim(*ylim)
		max_y = ylim[1]
		plt.plot(ylim, ylim, 'r-')
		plt.xlabel('old flux [mJy]')
		plt.ylabel('new flux [mJy]')
		name = ' '.join(cat.split('/')[-1].split('_')[:2])
		plt.title(name)
		outpath = cat.replace('combined_hdr_catalog.txt','xsc_new_vs_old_phot.png')
		plt.savefig(outpath, dpi=200)
		plt.close()
		print('\ncreated file: {}'.format(outpath))

	outfile = 'xsc_replaced.json'
	json.dump(replaced, open(outfile,'w'))
	print("\ncreated file: {}".format(outfile))
	print("\nremoved / added")
	for k,v in replaced.iteritems():
		print k.split('/')[-1], v['old'], v['new']
	m = np.mean([i['old']/float(i['new']) for i in replaced.values()])
	print("average ratio: {}".format(m))
	print("\nK mag and r_ext of sources with NaN photometry:")
	for i in find_files(bcdphot_out_path, "*xsc_nan_phot.csv"):
		reg = i.split('/')[-1]
		rec = np.recfromcsv(i)
		bad_id = rec.designation.tolist()
		bad_k = rec.k_m_k20fe.tolist()
		bad_r_ext = rec.r_ext.tolist()
		print reg
		print ("\tid\t\t\tKmag\tr_ext")
		if type(bad_id) is list:
			seq = sorted(zip(bad_id, bad_k, bad_r_ext), key=lambda x: x[0])
			for j,k,l in seq:
				print("\t{}\t{}\t{}".format(j,k,l))
		else:
			print("\t{}\t{}\t{}".format(bad_id, bad_k, bad_r_ext))


if __name__ == '__main__':
	bcdphot_out_path = '/Users/jlivings/data/spitzer_cirb_data/ace7/bcdphot_out_yesmask_nocentroid'
	mosaic_path = '/Users/jlivings/data/spitzer_cirb_data/mopex'
	run_xsc_phot(bcdphot_out_path, mosaic_path)