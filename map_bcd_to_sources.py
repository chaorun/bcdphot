#/usr/bin/env python

import os
import sys
import simplejson as json
import itertools
import glob
import multiprocessing

def map_to_sources(filepath):
	"""
	reads JSON bcd list file, gets the set of BCD files,
	then associates each with a set of RA/Dec coordinates
	"""
	sources = json.load(open(filepath))

	list2d = [i['files'] for i in sources]
	merged = list(itertools.chain.from_iterable(list2d))
	bcd_list = list(set(merged))
	bcd_list.sort()

	d = {}
	for i in bcd_list:
		d[i] = []
		for s in sources:
			if i in s['files']:
				d[i].append((s['ra'],s['dec']))

	outfilepath = filepath.replace('bcd_list.json','source_list.json')
	with open(outfilepath,'w') as w:
		json.dump(d,w,indent=4*' ')
	print('created file: '+outfilepath)

if __name__ == "__main__":
	filepath = sys.argv[1]
	map_to_sources(filepath)