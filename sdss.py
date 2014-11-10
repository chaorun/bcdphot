import os
import sys
import numpy as np
import pandas as pd
import urllib, urllib2
from prettyplotlib import plt
from astropy.io.votable import parse_single_table
from util import find_files
from util import spherematch
from util import great_circle_distance
from util import spz_jy_to_mags


def get_url(ra, dec, radius):

	"""
	Returns a URL for a HTTP query of the SDSS-DR9 catalog via ViZieR,
	centered on <ra>, <dec> (in decimal degrees), for a circle with
	with radius <radius> in degrees. URL will return data in XML/VOTable
	format.
	"""

	url = 'http://vizier.u-strasbg.fr/viz-bin/votable'
	p = {}
	p['-source'] = 'V/139'
	if np.sign(dec) > 0:
		p['-c'] = '{0:0.4f}+{1:0.4f}'.format(ra, dec)
	else:
		p['-c'] = '{0:0.4f}{1:0.4f}'.format(ra, dec)
	# s = np.sqrt(size)
	# p['-c.bd'] = '{0:.4f}x{1:.4f}'.format(s, s)
	p['-c.rd'] = '{}'.format(radius)
	# p['-out'] = '**'
	p['-out'] = ' '.join(['objID', 'RAJ2000', 'DEJ2000', 'cl'])
	p['-out.max'] = 999999999
	query = urllib.urlencode(p)
	query_url = url + "?" + query
	return query_url


def match_sdss(cat_path):
	for catfile in find_files(cat_path, "*merged.txt"):

		# read pipeline catalog
		print("\nreading catalog: {}".format(catfile))
		cat = pd.read_table(catfile, sep=' ')

		# retrieve SDSS data from ViZieR if not already downloaded
		ch = catfile.split('/')[-1].split('_')[1]
		outpath = catfile.replace('{}_merged.txt'.format(ch), 'sdss.vot')
		if not os.path.isfile(outpath):
			cntr_ra = np.median(cat.ra)
			cntr_dec = np.median(cat.dec)
			# get source from one corner of the mosaic to calculate radius
			c1 = (cat.ra.min(), cat.dec[cat.ra==cat.ra.min()].values[0])
			# make radius 10% bigger just to be on safe side
			radius = great_circle_distance(cntr_ra, cntr_dec, *c1) * 1.1
			url = get_url(cntr_ra, cntr_dec, radius)
			print("retrieving URL: {}".format(url))
			handler = urllib2.urlopen(url)
			raw = handler.read()
			with open(outpath,'wb') as f:
				f.write(raw)
			print("created file: {}".format(outpath))

		# parse VOTable
		print("reading VOTable: {}".format(outpath))
		table = parse_single_table(outpath)

		# if this is one of the southern hemisphere regions, delete and continue
		if table.array.size == 0:
			os.remove(outpath)
			print("outside of SDSS coverage")
			continue

		# make sure no missing data
		for name in table.array.dtype.names:
			assert table.array[name].mask.sum() == 0

		# get unmasked array
		sdss = table.array.data

		# make sure sky coverage is big enough
		assert sdss['RAJ2000'].min() < cat.ra.min()
		assert sdss['RAJ2000'].max() > cat.ra.max()
		assert sdss['DEJ2000'].min() < cat.dec.min()
		assert sdss['DEJ2000'].max() > cat.dec.max()

		# match to catalog
		assert cat.shape[0] < sdss.shape[0]
		tol = 2/3600.
		idx1, idx2, ds = spherematch(cat.ra, cat.dec, 
			sdss['RAJ2000'], sdss['DEJ2000'], tolerance = tol)
		print("matched {} out of {} sources with {} arcsec tolerance".format(ds.size, 
			cat.shape[0], tol*3600))

		# create vector of star/galaxy class (0=missing, 3=galaxy, 6=star)
		cl = np.zeros(cat.shape[0]).astype('int')
		cl[idx1] = sdss['cl'][idx2]

		# add the column to the dataset
		cat['cl'] = cl

		# write to new file
		outpath = catfile.replace('merged.txt', 'merged+sdss.txt')
		# fmt = ['%i']+['%0.8f']*2+['%.4e']*2+['%i']*2
		# hdr = ' '.join(names)+' cl'
		# np.savetxt(outpath, df.to_records(index=False), fmt=fmt, header=hdr)
		cat.to_csv(outpath, index=False, sep=' ', float_format='%.8f')
		print("created file: {}".format(outpath))


def plot_sdss(cat_path):
	for catfile in find_files(cat_path, "*merged+sdss.txt"):

		# for now ignore the channel 2 files
		if catfile.split('/')[-1].split('_')[1] != '1':
			continue

		print("\nreading catalog: {}".format(catfile))
		df = pd.read_table(catfile, sep=' ')

		# get rid of negative flux sources, if any
		df = df[df.flux > 0]

		# convert to magnitudes
		mags = spz_jy_to_mags(df.flux*1e-3, 1)

		# print counts per magnitude bin
		for i in range(10,15):
			sc = ((df.cl == 3) & (mags > i) & (mags < i+1)).sum()
			xc = ((df.xsc == 1) & (mags > i) & (mags < i+1)).sum() 
			msg = "{}th to {}th mag: {} SDSS galaxy sources, {} 2MASS XSC sources"
			print(msg.format(i, i+1, sc, xc))

		# print number of sources agreed upon
		agree = ((df.xsc == 1) & (df.cl == 3)).sum()
		disagree = ((df.xsc == 1) & (df.cl == 6)).sum()
		na = ((df.xsc == 1) & (df.cl == 0)).sum()
		msg = "{} 2MASS XSC sources classified as galaxies by SDSS"
		print(msg.format(agree))
		msg = "{} 2MASS XSC sources classified as stars by SDSS"
		print(msg.format(disagree))
		msg = "{} 2MASS XSC sources not matched to SDSS"
		print(msg.format(na))

		# plot normed histograms of 2MASS XSC and SDSS galaxy magnitudes
		xsc_gals = (mags > 10) & (mags < 15) & (df.xsc == 1)
		sdss_gals = (mags > 10) & (mags < 15) & (df.cl == 3)
		# mags[xsc_gals].hist(label='2MASS XSC', normed=True)
		# mags[sdss_gals].hist(label='SDSS galaxies', normed=True)
		plt.hist([mags[xsc_gals].values, mags[sdss_gals].values],
			bins=5, label=['2MASS', 'SDSS'])
		plt.xlabel('IRAC1 [mag]')
		plt.ylabel('Number Count')
		reg = catfile.split('/')[-1].split('_')[0]
		plt.title('{} Extended Sources / Galaxies'.format(reg))
		plt.legend(loc=2)
		name = '{}_2mass_xsc_vs_sdss_hist.png'.format(reg)
		outpath = '/'.join(catfile.split('/')[:-1]+[name])
		plt.savefig(outpath, dpi=100)
		plt.close()
		print("created file: {}".format(outpath))


if __name__ == '__main__':
	cat_path = '/Users/jlivings/data/spitzer_cirb_data/subarray/finalphot10_linefits_catalogs'
	match_sdss(cat_path)
	plot_sdss(cat_path)