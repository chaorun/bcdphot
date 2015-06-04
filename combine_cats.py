import os
import sys
import numpy as np
import pandas as pd
bcdphot_dir = '/Users/jlivings/data/spikes/bcdphot'
sys.path.append(bcdphot_dir)
from util import spherematch, spz_jy_to_mags, get_complement
import matplotlib.pyplot as plt


base_dir = '/Users/jlivings/data/spikes/railspikes_database_prep'
cat_dir = os.path.join(base_dir, 'complete_cats')

# step 1: combine individual tile catalogs
# ------------------------------------------

filepaths = [os.path.join(cat_dir, f) for f in os.listdir(cat_dir)]
filepaths = filter(lambda x: '.csv' in x, filepaths)
filepaths = filter(lambda x: 'epoch1' not in x, filepaths)
df1 = pd.read_csv(filepaths[0])
ra1 = df1[['ra_1', 'ra_2']].mean(axis=1)
dec1 = df1[['dec_1', 'dec_2']].mean(axis=1)

matched = 0
matched_ra, matched_dec = np.array([]), np.array([])
for fp in filepaths[1:]:

	df2 = pd.read_csv(fp)
	ra2 = df2[['ra_1', 'ra_2']].mean(axis=1)
	dec2 = df2[['dec_1', 'dec_2']].mean(axis=1)

	if ra1.size < ra2.size:
		idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2)
	else:
		idx2, idx1, ds = spherematch(ra2, dec2, ra1, dec1)

	# keep track of how many got matched and what their (mean) coordinates are
	matched += idx1.size
	new_ra = np.c_[df1[['ra_1', 'ra_2']].loc[idx1].values, df2[['ra_1', 'ra_2']].loc[idx2].values].mean(axis=1)
	new_dec = np.c_[df1[['dec_1', 'dec_2']].loc[idx1].values, df2[['dec_1', 'dec_2']].loc[idx2].values].mean(axis=1)
	matched_ra = np.append(matched_ra, new_ra)
	matched_dec = np.append(matched_dec, new_dec)
	# combine the matched sources into single sources before adding
	df1_subset = df1.loc[idx1]
	df2_subset = df2.loc[idx2]
	df1_subset.index = range(df1_subset.shape[0])
	df2_subset.index = range(df2_subset.shape[0])
	concat = pd.concat([df1_subset, df2_subset])
	concat_mean = concat.groupby(concat.index).mean()
	# for now just quick and dirty calculation of uncertainties -- completely
	# correct implementation will require making use of systematic uncertainty values,
	# i.e. unc_phot ** 2 = unc_tot ** 2 - unc_sys ** 2
	concat_mean.unc_1 /= np.sqrt(2)
	concat_mean.unc_2 /= np.sqrt(2)
	concat_sum = concat.groupby(concat.index).sum()
	concat_mean.n_obs_1 = concat_sum.n_obs_1
	concat_mean.n_obs_2 = concat_sum.n_obs_2
	# update df1[idx1] with these new values
	not_idx1 = get_complement(idx1, ra1.size)
	# data = np.r_[df1.loc[not_idx1].values, concat_mean.values]
	df1 = pd.concat([df1.loc[not_idx1], concat_mean])
	df1.index = range(df1.shape[0])

	# add new rows that don't have a match by position
	not_idx2 = get_complement(idx2, ra2.size)
	df1 = pd.concat([df1, df2.loc[not_idx2]])
	df1.index = range(df1.shape[0])

	# update ra1 and dec1 vectors for next iteration
	ra1 = df1[['ra_1', 'ra_2']].mean(axis=1)
	dec1 = df1[['dec_1', 'dec_2']].mean(axis=1)

print("{} matched sources combined into single measurements".format(matched))

# prepare output dataset
ra = df1[['ra_1', 'ra_2']].mean(axis=1)
dec = df1[['dec_1', 'dec_2']].mean(axis=1)
i1 = spz_jy_to_mags(df1['flux_1'] * 1e-3, 1)
i1unc = 1.08 * df1['unc_1'] / df1['flux_1']
i2 = spz_jy_to_mags(df1['flux_2'] * 1e-3, 2)
i2unc = 1.08 * df1['unc_2'] / df1['flux_2']
df = pd.DataFrame(dict(ra=ra, dec=dec, i1_mag=i1, i1_unc=i1unc, 
	i2_mag=i2, i2_unc=i2unc, i1_obs=df1['n_obs_1'], i2_obs=df1['n_obs_2']))
df = df[['ra', 'dec', 'i1_mag', 'i1_unc', 'i2_mag', 'i2_unc', 'i1_obs', 'i2_obs']]
outpath	= os.path.join(base_dir, 'combined_cats.csv')
df.to_csv(outpath, index=False)

# make plots
plt.hist2d(ra, dec, bins=180, cmap=plt.cm.hot)
outpath = os.path.join(base_dir, 'combined_cats.pdf')
plt.savefig(outpath, dpi=120)
plt.close()
plt.hist2d(matched_ra, matched_dec, bins=180, cmap=plt.cm.hot)
xlim, ylim = plt.xlim(), plt.ylim()
outpath = os.path.join(base_dir, 'tile_overlap_matched_sources.pdf')
plt.savefig(outpath, dpi=120)
plt.close()
plt.hist2d(ra, dec, bins=180, cmap=plt.cm.hot)
plt.xlim(xlim)
plt.ylim(ylim)
outpath = os.path.join(base_dir, 'combined_cats_zoom.pdf')
plt.savefig(outpath, dpi=120)
plt.close()

# step 2: match to the kic koi ukirtj data
# ------------------------------------------
# outpath	= os.path.join(base_dir, 'combined_cats.csv')
# df = pd.read_csv(outpath)
kic = pd.read_csv('/Users/jlivings/data/spikes/other_catalogs/kic_kois_ukirtj.csv')
ra1, dec1 = kic.kic_ra, kic.kic_dec
ra2, dec2 = df.ra, df.dec
if ra1.size < ra2.size:
	idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2, tolerance=2/3600.)
else:
	idx2, idx1, ds = spherematch(ra2, dec2, ra1, dec1, tolerance=2/3600.)
kic_matched = kic.loc[idx1]
df_matched = df.loc[idx2]
for i in df.columns[2:]:
	kic_matched[i] = df_matched[i].values
outpath = os.path.join(base_dir, 'combined_cats+kic_kois_ukirtj.csv')
kic_matched.to_csv(outpath, index=False)
outpath	= os.path.join(base_dir, 'combined_cats_matched_to_kois.pdf')
plt.hexbin(kic_matched.kic_ra, kic_matched.kic_dec)
plt.savefig(outpath, dpi=120)
plt.close()

# now add two more columns for colors
kic_matched['i1i2color'] = kic_matched['i1_mag'] - kic_matched['i2_mag']
kic_matched['ji2color'] = kic_matched['ukirt_j'] - kic_matched['i2_mag']
outpath = os.path.join(base_dir, 'combined_cats+kic_kois_ukirtj_color.csv')
kic_matched.to_csv(outpath, index=False)


# step 3: match to the kic kob ukirtj data
# ------------------------------------------

kic = pd.read_csv('/Users/jlivings/data/spikes/other_catalogs/kic_kobs_ukirtj.csv')
ra1, dec1 = kic.kic_ra, kic.kic_dec
ra2, dec2 = df.ra, df.dec
if ra1.size < ra2.size:
	idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2, tolerance=2/3600.)
else:
	idx2, idx1, ds = spherematch(ra2, dec2, ra1, dec1, tolerance=2/3600.)
kic_matched = kic.loc[idx1]
df_matched = df.loc[idx2]
for i in df.columns[2:]:
	kic_matched[i] = df_matched[i].values
outpath = os.path.join(base_dir, 'combined_cats+kic_kobs_ukirtj.csv')
kic_matched.to_csv(outpath, index=False)
outpath	= os.path.join(base_dir, 'combined_cats_matched_to_kobs.pdf')
plt.hexbin(kic_matched.kic_ra, kic_matched.kic_dec)
plt.savefig(outpath, dpi=120)
plt.close()

# now add two more columns for colors
kic_matched['i1i2color'] = kic_matched['i1_mag'] - kic_matched['i2_mag']
kic_matched['ji2color'] = kic_matched['ukirt_j'] - kic_matched['i2_mag']
outpath = os.path.join(base_dir, 'combined_cats+kic_kobs_ukirtj_color.csv')
kic_matched.to_csv(outpath, index=False)


# just plots

def plot_simple_pos(ra, dec, outpath):
	# plt.hexbin(kois.kic_ra, kois.kic_dec, cmap=plt.cm.Greys)
	fig = plt.figure(frameon=False)
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	ax.set_axis_off()
	fig.add_axes(ax)
	ax.plot(ra, dec, 'k.', ms=1)
	plt.xlim(ra.min(), ra.max())
	plt.ylim(dec.min(), dec.max())
	plt.savefig(outpath, dpi=120)
	plt.close()

kois = pd.read_csv('combined_cats+kic_kois_ukirtj.csv')
kobs = pd.read_csv('combined_cats+kic_kobs_ukirtj.csv')

outpath	= os.path.join(base_dir, 'KOIs.png')
plot_simple_pos(kois.kic_ra, kois.kic_dec, outpath)

outpath	= os.path.join(base_dir, 'KOBs.png')
plot_simple_pos(kobs.kic_ra, kobs.kic_dec, outpath)
