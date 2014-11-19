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
	Returns a URL for a HTTP query of the WISE catalog via ViZieR,
	centered on <ra>, <dec> (in decimal degrees), for a circle with
	with radius <radius> in degrees. URL will return data in XML/VOTable
	format.
	"""

	url = 'http://vizier.u-strasbg.fr/viz-bin/votable'
	p = {}
	p['-source'] = 'II/311'
	if np.sign(dec) > 0:
		p['-c'] = '{0:0.4f}+{1:0.4f}'.format(ra, dec)
	else:
		p['-c'] = '{0:0.4f}{1:0.4f}'.format(ra, dec)
	# s = np.sqrt(size)
	# p['-c.bd'] = '{0:.4f}x{1:.4f}'.format(s, s)
	p['-c.rd'] = '{}'.format(radius)
	# p['-out'] = '**'
	p['-out'] = ' '.join(['WISE', 'RAJ2000', 'DEJ2000', 'W1mag', 
		'e_W1mag', 'W2mag', 'e_W2mag'])
	p['-out.max'] = 999999999
	query = urllib.urlencode(p)
	query_url = url + "?" + query
	return query_url


def match_wise(cat_path, sdss=True):
	if sdss:
		search_pattern = "*merged+sdss.txt"
	else:
		search_pattern = "*merged.txt"

	for catfile in find_files(cat_path, search_pattern):

		# read pipeline catalog
		print("\nreading catalog: {}".format(catfile))
		cat = pd.read_table(catfile, sep=' ')

		# retrieve WISE data from ViZieR if not already downloaded
		ch = catfile.split('/')[-1].split('_')[1]
		if sdss:
			outpath = catfile.replace('{}_merged+sdss.txt'.format(ch), 'wise.vot')
		else:
			outpath = catfile.replace('{}_merged.txt'.format(ch), 'wise.vot')
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
			print("no WISE coverage")
			continue

		# get unmasked array
		wise = table.array.data

		# make sure sky coverage is big enough
		assert wise['RAJ2000'].min() < cat.ra.min()
		assert wise['RAJ2000'].max() > cat.ra.max()
		assert wise['DEJ2000'].min() < cat.dec.min()
		assert wise['DEJ2000'].max() > cat.dec.max()

		# match to catalog
		tol = 2/3600.
		if cat.shape[0] < wise.shape[0]:
			idx1, idx2, ds = spherematch(cat.ra, cat.dec, 
				wise['RAJ2000'], wise['DEJ2000'], tolerance = tol)
		else:
			idx2, idx1, ds = spherematch(wise['RAJ2000'], wise['DEJ2000'],
				cat.ra, cat.dec, tolerance = tol)
		print("matched {} out of {} sources with {} arcsec tolerance".format(ds.size, 
			cat.shape[0], tol*3600))

		# add WISE to the catalog
		if ch == '1':
			cat['W1mag'] = np.repeat(np.nan, cat.shape[0])
			cat['e_W1mag'] = np.repeat(np.nan, cat.shape[0])
			cat['W1mag'][idx1] = wise['W1mag'][idx2]
			cat['e_W1mag'][idx1] = wise['e_W1mag'][idx2]
		elif ch == '2':
			cat['W2mag'] = np.repeat(np.nan, cat.shape[0])
			cat['e_W2mag'] = np.repeat(np.nan, cat.shape[0])
			cat['W2mag'][idx1] = wise['W2mag'][idx2]
			cat['e_W2mag'][idx1] = wise['e_W2mag'][idx2]
		else:
			print("unexpected error adding WISE data")

		# write to new file
		outpath = catfile.replace('.txt', '+wise.csv')
		# fmt = ['%i']+['%0.8f']*2+['%.4e']*2+['%i']*2
		# hdr = ' '.join(names)+' cl'
		# np.savetxt(outpath, df.to_records(index=False), fmt=fmt, header=hdr)
		cat.to_csv(outpath, index=False, float_format='%.8f')
		print("created file: {}".format(outpath))

def plot_wise(cat_path):

	for catfile in find_files(cat_path, "*merged+wise.csv"):

		print("\nreading catalog: {}".format(catfile))
		df = pd.read_csv(catfile)

		# convert to magnitudes
		nbadflux = (df.flux <= 0).sum()
		try:
			assert nbadflux == 0
		except:
			print("warning: {} negative flux source(s)".format(nbadflux))
		ch = catfile.split('/')[-1].split('_')[1]
		mags = spz_jy_to_mags(df.flux*1e-3, float(ch))
		if ch == '1':
			plt.scatter(df.W1mag, mags)
			plt.xlabel('W1 [mag]')
			plt.ylabel('I1 [mag]')
		elif ch == '2':
			plt.scatter(df.W2mag, mags)
			plt.xlabel('W2 [mag]')
			plt.ylabel('I2 [mag]')
		ax = plt.gca()
		xlim, ylim = ax.get_xlim(), ax.get_ylim()
		plt.plot([-5, ylim[1]*2], [-5, ylim[1]*2], 'r-')
		ax.set_xlim(xlim) ; ax.set_ylim(ylim)
		reg = catfile.split('/')[-1].split('_')[0]
		name = '{}_{}_IRAC_vs_WISE.png'.format(reg, ch)
		outpath = '/'.join(catfile.split('/')[:-1]+[name])
		plt.savefig(outpath, dpi=120)
		plt.close()


def match_cats(df1, df2, tol=2/3600.):
	assert 'ra' in df1.columns
	assert 'dec' in df1.columns
	assert 'ra' in df2.columns
	assert 'dec' in df2.columns
	if df1.shape[0] < df2.shape[0]:
		idx1, idx2, ds = spherematch(df1.ra, df1.dec, 
			df2.ra, df2.dec, tolerance=tol)
	else:
		idx2, idx1, ds = spherematch(df2.ra, df2.dec, 
			df1.ra, df1.dec, tolerance=tol)
	return idx1, idx2


def plot(x, y, outpath, xlabel, ylabel, plot_style, plot_type):
	if plot_type == 'mag-mag':
		xlim = (10, 18)
		ylim = (10, 18)
	elif plot_type == 'color-mag':
		xlim = (10, 18)
		ylim = (-1, 1)
	elif plot_type == 'color-color':
		xlim = (-1, 1)
		ylim = (-1, 1)
	else:
		raise(ValueError("plot_type should be one of ['mag-mag', 'color-mag', 'color-color'] "))
	isinrange = lambda a,b: (a>=b[0]) & (a<=b[1])
	g = isinrange(x, xlim) & isinrange(y, ylim)
	if plot_style == 'scatter':
		plt.scatter(x[g], y[g])
	elif plot_style == 'hexbin':
		plt.hexbin(x[g], y[g])
	elif plot_style == 'hist2d':
		plt.hist2d(x[g], y[g], bins=100)
	else:
		raise(ValueError("plot_style should be one of ['scatter', 'hexbin', 'hist2d'] "))
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	ax = plt.gca()
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	plt.savefig(outpath, dpi=120)
	plt.close()
	print("created file: {}".format(outpath))



def plot_spz_vs_wise(cat_path, plot_style='scatter'):

	ch1 = list(find_files(cat_path, "*merged+wise.csv"))[::2]
	ch2 = list(find_files(cat_path, "*merged+wise.csv"))[1::2]

	for ch1, ch2 in zip(ch1, ch2):

		reg1 = ch1.split('/')[-1].split('_')[0]
		reg2 = ch2.split('/')[-1].split('_')[0]
		assert reg1 == reg2

		print("\nreading catalog: {}".format(ch1))
		print("reading catalog: {}".format(ch2))
		df1 = pd.read_csv(ch1)
		df2 = pd.read_csv(ch2)

		# convert to magnitudes
		mags1 = spz_jy_to_mags(df1.flux*1e-3, 1)
		mags2 = spz_jy_to_mags(df2.flux*1e-3, 2)

		# match ch1 / ch2
		idx1, idx2 = match_cats(df1, df2, tol=2/3600.)

		# save matched catalogs
		matched1 = df1.loc[idx1]
		matched2 = df2.loc[idx2]
		ch1_cols = [i+'_1' for i in df1.columns.tolist()]
		ch2_cols = [i+'_2' for i in df2.columns.tolist()]
		matched1.columns = ch1_cols	
		matched2.columns = ch2_cols
		# matched = pd.concat([matched1, matched2], 1, ignore_index=True)	# weird error
		matched = np.concatenate([matched1.values, matched2.values], 1)
		df_matched = pd.DataFrame(matched, columns=ch1_cols+ch2_cols)
		df_matched['I1'] = mags1[idx1].values
		df_matched['I2'] = mags2[idx2].values
		outpath = '/'.join(ch1.split('/')[:-1])+'/{}_2ch_matched.csv'.format(reg1)
		df_matched.to_csv(outpath, index=False, float_format='%.8f')
		print("created file: {}".format(outpath))

		# plot I1-I2 vs. W1-W2
		color1 = df1.W1mag[idx1].values - df2.W2mag[idx2].values
		color2 = mags1[idx1].values - mags2[idx2].values
		name = '{}_I1-I2_vs_W1-W2_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(color1, color2, outpath, 'W1-W2 [mag]', 'I1-I2 [mag]', 
			plot_style=plot_style, plot_type='color-color')

		# plot I1-W1 vs. I2-W2
		color1 = mags1[idx1].values - df1.W1mag[idx1].values
		color2 = mags2[idx2].values - df2.W2mag[idx2].values
		name = '{}_I1-W1_vs_I2-W2_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(color1, color2, outpath, 'I1-W1 [mag]', 'I2-W2 [mag]', 
			plot_style=plot_style, plot_type='color-color')

		# plot spz color-magnitude diagrams
		color = mags1[idx1].values - mags2[idx2].values
		mags = mags1[idx1].values
		name = '{}_I1_vs_I1-I2_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(mags, color, outpath, 'I1 [mag]', 'I1-I2 [mag]', 
			plot_style=plot_style, plot_type='color-mag')

		# plot wise color-magnitude diagrams
		color = df1.W1mag[idx1].values - df2.W2mag[idx2].values
		mags = df1.W1mag[idx1].values
		name = '{}_W1_vs_W1-W2_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(mags, color, outpath, 'W1 [mag]', 'W1-W2 [mag]', 
			plot_style=plot_style, plot_type='color-mag')
	
		# plot I1 vs I2
		mags1_matched = mags1[idx1].values
		mags2_matched = mags2[idx2].values
		name = '{}_I1_vs_I2_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(mags1_matched, mags2_matched, outpath, 'I1 [mag]', 'I2 [mag]', 
			plot_style=plot_style, plot_type='mag-mag')


def plot_spz_vs_wise_sdss_class(cat_path, plot_style='scatter'):

	ch1 = list(find_files(cat_path, "*merged+sdss+wise.csv"))[::2]
	ch2 = list(find_files(cat_path, "*merged+sdss+wise.csv"))[1::2]

	for ch1, ch2 in zip(ch1, ch2):

		reg1 = ch1.split('/')[-1].split('_')[0]
		reg2 = ch2.split('/')[-1].split('_')[0]
		assert reg1 == reg2

		print("\nreading catalog: {}".format(ch1))
		print("reading catalog: {}".format(ch2))
		df1 = pd.read_csv(ch1)
		df2 = pd.read_csv(ch2)

		# convert to magnitudes
		mags1 = spz_jy_to_mags(df1.flux*1e-3, 1)
		mags2 = spz_jy_to_mags(df2.flux*1e-3, 2)

		# match ch1 / ch2
		idx1, idx2 = match_cats(df1, df2, tol=2/3600.)

		# save matched catalogs
		matched1 = df1.loc[idx1]
		matched2 = df2.loc[idx2]
		ch1_cols = [i+'_1' for i in df1.columns.tolist()]
		ch2_cols = [i+'_2' for i in df2.columns.tolist()]
		matched1.columns = ch1_cols	
		matched2.columns = ch2_cols
		# matched = pd.concat([matched1, matched2], 1, ignore_index=True)	# weird error
		matched = np.concatenate([matched1.values, matched2.values], 1)
		df_matched = pd.DataFrame(matched, columns=ch1_cols+ch2_cols)
		df_matched['I1'] = mags1[idx1].values
		df_matched['I2'] = mags2[idx2].values
		outpath = '/'.join(ch1.split('/')[:-1])+'/{}_2ch_matched+sdss.csv'.format(reg1)
		df_matched.to_csv(outpath, index=False, float_format='%.8f')
		print("created file: {}".format(outpath))

		# identify SDSS galaxies and stars
		galaxies = (df1.cl[idx1].values == 3) & (df2.cl[idx2].values == 3)
		stars = (df1.cl[idx1].values == 6) & (df2.cl[idx2].values == 6)

		# plot I1-I2 vs. W1-W2
		color1 = df1.W1mag[idx1].values - df2.W2mag[idx2].values
		color2 = mags1[idx1].values - mags2[idx2].values
		# galaxies
		name = '{}_I1-I2_vs_W1-W2_galaxies_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(color1[galaxies], color2[galaxies], outpath, 'W1-W2 [mag]', 'I1-I2 [mag]', 
			plot_style=plot_style, plot_type='color-color')
		# stars
		outpath = '/'.join(ch1.split('/')[:-1]+[name]).replace('galaxies', 'stars')
		plot(color1[stars], color2[stars], outpath, 'W1-W2 [mag]', 'I1-I2 [mag]', 
			plot_style=plot_style, plot_type='color-color')

		# plot I1-W1 vs. I2-W2
		color1 = mags1[idx1].values - df1.W1mag[idx1].values
		color2 = mags2[idx2].values - df2.W2mag[idx2].values
		# galaxies
		name = '{}_I1-W1_vs_I2-W2_galaxies_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(color1[galaxies], color2[galaxies], outpath, 'I1-W1 [mag]', 'I2-W2 [mag]', 
			plot_style=plot_style, plot_type='color-color')
		# stars
		outpath = '/'.join(ch1.split('/')[:-1]+[name]).replace('galaxies', 'stars')
		plot(color1[stars], color2[stars], outpath, 'I1-W1 [mag]', 'I2-W2 [mag]', 
			plot_style=plot_style, plot_type='color-color')

		# plot spz color-magnitude diagrams
		color = mags1[idx1].values - mags2[idx2].values
		mags = mags1[idx1].values
		# galaxies
		name = '{}_I1_vs_I1-I2_galaxies_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(mags[galaxies], color[galaxies], outpath, 'I1 [mag]', 'I1-I2 [mag]', 
			plot_style=plot_style, plot_type='color-mag')
		# stars
		outpath = '/'.join(ch1.split('/')[:-1]+[name]).replace('galaxies', 'stars')
		plot(mags[stars], color[stars], outpath, 'I1 [mag]', 'I1-I2 [mag]', 
			plot_style=plot_style, plot_type='color-mag')

		# plot wise color-magnitude diagrams
		color = df1.W1mag[idx1].values - df2.W2mag[idx2].values
		mags = df1.W1mag[idx1].values
		# galaxies
		name = '{}_W1_vs_W1-W2_galaxies_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(mags[galaxies], color[galaxies], outpath, 'W1 [mag]', 'W1-W2 [mag]', 
			plot_style=plot_style, plot_type='color-mag')
		# stars
		outpath = '/'.join(ch1.split('/')[:-1]+[name]).replace('galaxies', 'stars')
		plot(mags[stars], color[stars], outpath, 'W1 [mag]', 'W1-W2 [mag]', 
			plot_style=plot_style, plot_type='color-mag')
	
		# plot I1 vs I2
		mags1_matched = mags1[idx1].values
		mags2_matched = mags2[idx2].values
		# galaxies
		name = '{}_I1_vs_I2_galaxies_plot_style.png'.format(reg1)
		name = name.replace('plot_style', plot_style)
		outpath = '/'.join(ch1.split('/')[:-1]+[name])
		plot(mags1_matched[galaxies], mags2_matched[galaxies], outpath, 'I1 [mag]', 'I2 [mag]', 
			plot_style=plot_style, plot_type='mag-mag')
		# stars
		outpath = '/'.join(ch1.split('/')[:-1]+[name]).replace('galaxies', 'stars')
		plot(mags1_matched[stars], mags2_matched[stars], outpath, 'I1 [mag]', 'I2 [mag]', 
			plot_style=plot_style, plot_type='mag-mag')


if __name__ == '__main__':
	cat_path = '/Users/jlivings/data/spitzer_cirb_data/subarray/finalphot10_linefits_catalogs'
	match_wise(cat_path, sdss=True)
	plot_wise(cat_path)
	plot_spz_vs_wise(cat_path, plot_style='scatter')
	plot_spz_vs_wise(cat_path, plot_style='hexbin')
	plot_spz_vs_wise(cat_path, plot_style='hist2d')
	plot_spz_vs_wise_sdss_class(cat_path, plot_style='scatter')
	plot_spz_vs_wise_sdss_class(cat_path, plot_style='hexbin')
	plot_spz_vs_wise_sdss_class(cat_path, plot_style='hist2d')