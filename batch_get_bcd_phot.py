import multiprocessing
import glob
from get_bcd_phot import *

source_list_paths = glob.glob('bcd_dirs/*/*/*source_list.json')

def process(source_list_path):
	data_dir = source_list_path.split('source_list.json')[0]
	cbuncpaths = glob.glob('unzipdirs/*/*/*/*cbunc.fits')
	gross_lst = get_gross_list(source_list_path,cbuncpaths)
	gross_arr = np.array(gross_lst).astype(np.float)
	np.savetxt(data_dir+'gross_arr.txt',gross_arr)
	print('created file: '+data_dir+'gross_arr.txt')
	phot_groups = get_phot_groups(gross_arr,data_dir)

ncpus = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=ncpus)
print "using %i CPUs" % ncpus

pool.map(process,source_list_paths)