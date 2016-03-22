import os, sys
import numpy as np
import pandas as pd
from util import spherematch
from util import find_files
from util import ordinary_least_squares
import matplotlib.pyplot as plt
import statsmodels.api as sm

# work dir: ~/data/spitzer_cirb_data/subarray

def cleanup_sub(vg_dir):
	"""
	clean up the photometry data from Varoujan
	"""
	new_dir = vg_dir+'_clean'
	os.mkdir(new_dir)
	phot_vg_files = filter(lambda x: '.txt' in x, os.listdir(vg_dir))
	# phot_vg_phottot_files = filter(lambda x: 'phottot' in x, phot_vg_files)
	for f in phot_vg_files:
		df = pd.read_table(vg_dir+'/'+f,
			names = ['id','ra','dec','flux','unc','x','y','flux_uncor'],
			delim_whitespace=True)
		starnums, dithers = zip(*[i.split('_')[1:4:2] for i in df.id])
		df['id'] = [int(i) for i in starnums]
		df['dither'] = [int(i) for i in dithers]
		sorted_df = df.sort(['id','dither'])
		# new: remove the aperture correction applied by varoujan to the uncertainties
		ch = f.split('-')[2]
		if ch == 'ch1':
			sorted_df['unc'] /= 1.205
		elif ch == 'ch2':
			sorted_df['unc'] /= 1.221
		else:
			raise(TypeError("unexpected channel"))
		fnew = '_'.join(f.split('-')[::2])+'_raw.csv'
		sorted_df.to_csv(new_dir+'/'+fnew, index=False, float_format='%.8f')
		# also calculate mean RA/Dec, flux, and quadrature sum uncertainty
		grouped = sorted_df.groupby('id')
		agg = grouped[['ra','dec','flux']].aggregate(np.median)
		quadsum = grouped['unc'].aggregate(lambda x: np.sqrt(np.sum(x**2)))
		agg['unc'] = quadsum
		fnew = '_'.join(f.split('-')[::2])+'_agg.csv'
		agg.to_csv(new_dir+'/'+fnew, index=True, float_format='%.8f')

def get_cutoff(channel):
	if float(channel) == 1:
		cutoff = 250
	elif float(channel) == 2:
		cutoff = 350
	else:
		raise ValueError('unexpected channel value')
	return cutoff

def fit_line(df, channel):
	cutoff = get_cutoff(channel)
	idx = (df['sub_flux'] > 0) & (df['sub_flux'] < cutoff)
	X = df['hdr_flux'][idx]
	y = df['sub_flux'][idx]
	# model = sm.OLS(y,X)
	# result = model.fit()
	# slope = result.params[0]
	slope = ordinary_least_squares(y, X)
	return slope

def merge_subarray(vg_dir, bcdphot_dir):
	out_dir = vg_dir.replace('clean','plots_catalogs')
	os.mkdir(out_dir)

	hdr_files = find_files(bcdphot_dir, '*combined_hdr_*xsc_cor.txt')
	# hdr_file = list(hdr_files)[0]
	for hdr_file in hdr_files:
		reg, ch = hdr_file.split('/')[-1].split('_')[:2]
		sub_file = '/'.join([vg_dir, "d{}_ch{}_agg.csv".format(reg, ch)])

		hdr_names = open(hdr_file).readline().split()[1:]
		hdr = np.recfromtxt(hdr_file, names=hdr_names)
		sub = np.recfromcsv(sub_file)
		# sub.flux *= 1e-3	# convert from uJy to mJy

		idx1, idx2, ds = spherematch(sub.ra, sub.dec, hdr.ra, hdr.dec, tolerance=3/3600.)
		df = pd.DataFrame({'sub_flux': sub.flux[idx1], 'hdr_flux':hdr.flux[idx2]})

		slope = fit_line(df, int(ch))
		with open("{}/linefits.txt".format(out_dir),'a') as f:
			f.write("{} {} {}\n".format(reg, ch, slope))

		fig = df.plot(x='hdr_flux',y='sub_flux', kind='scatter')
		fig.plot([0, fig.get_xlim()[1]], [0, slope * fig.get_xlim()[1]], 'r-')
		fig.set_title("region {} channel {}".format(reg, ch))
		fig.text(fig.get_xlim()[1]*0.2, fig.get_ylim()[1]*0.8, 
			"slope: {0:3f}".format(slope), fontsize=24)
		plt.savefig("{}/{}_{}_linefit.png".format(out_dir, reg, ch), dpi=300)
		plt.close()

		# now save the (uncorrected) matched data to disk
		sub_matched = pd.DataFrame.from_records(sub[idx1])
		# rename the columns
		cols = sub_matched.columns.tolist()
		cols_new = ['sub_'+i for i in cols]
		sub_matched.columns = cols_new
		# set hdr_matched dataframe index equal to sub_matched index, this is
		# necessary for concatenation using pandas.concat
		hdr_matched = pd.DataFrame.from_records(hdr[idx2]).set_index(sub_matched.index)
		# rename the columns
		cols = hdr_matched.columns.tolist()
		cols_new = ['hdr_'+i for i in cols]
		hdr_matched.columns = cols_new
		# concatenate
		concat = pd.concat([ sub_matched, hdr_matched ], 1)
		# # convert subarray flux to mJy
		# concat.sub_flux = concat.sub_flux*1e3
		# concat.sub_unc = concat.sub_unc*1e3
		concat.to_csv("{}/{}_{}_hdr_vs_sub.csv".format(out_dir, reg, ch), 
			index=False, float_format='%.8f')

		# now correct all the subarray flux values with the slope
		sub.flux /= slope

		# now merge hdr and subarray into one dataset:
		# want to keep all the hdr photometry that is not saturated, and
		# keep only the subarray photometry above the hdr saturation limit
		cutoff = get_cutoff(ch)
		bad = hdr.flux > cutoff
		hdr_subset = pd.DataFrame.from_records(hdr[~bad])
		bad = sub.flux < cutoff
		sub_subset = pd.DataFrame.from_records(sub[~bad])
		# add n_obs column to subarray data so it has same format as hdr
		sub_subset['n_obs'] = 4
		# add column indicating whether if it came from subarray
		hdr_subset['sub'] = np.zeros(hdr_subset.shape[0]).astype(int)
		sub_subset['sub'] = np.ones(sub_subset.shape[0]).astype(int)
		# concatenate them
		concat = pd.concat([ hdr_subset, sub_subset ], 0, ignore_index=True)
		# get rid of the 'id' field since it is no longer relevant
		# but add a column indicating if it was a 2MASS XSC measurement
		concat['xsc'] = np.zeros(concat.shape[0]).astype(int)
		concat.xsc[concat.id < 1] = 1
		concat = concat.drop('id', 1)
		# apply 1% flux reduction to correct for stray light (only to >100 mJy sources)
		concat.flux[concat.flux > 100] *= 0.99
		concat.unc[concat.flux > 100] *= 0.99
		# write to disk
		concat.to_csv("{}/{}_{}_merged.txt".format(out_dir, reg, ch), 
			index=False, sep=' ', float_format='%.8f')


if __name__ == '__main__':
	vg_dir = '/Users/jlivings/data/spitzer_cirb_data/subarray/finalphot10'
	cleanup_sub(vg_dir)

	vg_dir += '_clean'
	bcdphot_dir = '/Users/jlivings/data/spitzer_cirb_data/ace7/bcdphot_out_yesmask_nocentroid'
	merge_subarray(vg_dir, bcdphot_dir)