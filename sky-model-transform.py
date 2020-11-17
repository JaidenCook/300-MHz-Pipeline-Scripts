#!/usr/bin/python

__author__ = "Jaiden Cook"
__credits__ = ["Jaiden Cook"]
__version__ = "1.0.1"
__maintainer__ = "Jaiden Cook"
__email__ = "Jaiden.Cook@student.curtin.edu"

# Version 2 of the direct beam skymodel pipeline.

# Generic stuff:
import os
import sys
import time
from datetime import datetime
import glob
import shutil
import re
from math import pi

# Array stuff:
import numpy as np

# Plotting stuff:
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Multiprocessing stuff:
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

# Scipy stuff:
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.special import gamma, digamma
from scipy import stats
from scipy.stats import chi2
from scipy.stats import f

# Astropy stuff:
from astropy import wcs
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import Table,Column
from astropy.io.votable import writeto as writetoVO

sys.path.append(os.path.abspath("/home/jaiden/Documents/Masters_Project/bin"))
from JOOF import *

if __name__ == "__main__":
	# Timing the script duration.
	start0 = time.time()
	
	# Parser options:
	from optparse import OptionParser
	
	usage="Usage: %prog [options]\n"
	parser = OptionParser(usage=usage)
	parser.add_option('--catalogue',dest="catalogue",default=None,help="Input catalogue.")
	parser.add_option('--old_freq',dest="old_freq",default=None,help="Reference frequency of the old catalogue.")
	parser.add_option('--nu_freq',dest="nu_freq",default=None,help="Reference frequency for the ouptut catalogue.")
	parser.add_option('--nu_catalogue',dest="nu_catalogue",default=None,help="Name of new catalogue.")
	
	# Add option for primary beam model inputs -- make default a typical GLEAMY2 repository location
	(options, args) = parser.parse_args()

	try:
		catalogue = str(options.catalogue)
		old_freq = float(options.old_freq)
		nu_freq = float(options.nu_freq)
	except IOError:
		print "Incorrect I/O"
		system.exit(0)

	# Loading in the sky-model
	skymodel=fits.getdata("{0}".format(catalogue))
	
	# Turning the sky-model into an Astropy table.
	t=Table(skymodel)
	
	# Subsetting the old coefficients for testing purposes.
	t_co_old = np.array(t['coefficients'])
	
	# log_poly_co_map(a,x_a,x_b)
	t['coefficients'] = [(log_poly_co_map(t['coefficients'][i],old_freq,nu_freq)) for i in range(len(t))]
	
	# Renaming the flux density column to the appropriate name.
	t.rename_column('Fint{0}'.format(int(round(old_freq))),'Fint{0}'.format(int(round(nu_freq))))

	# Setting the new integrated flux densities.
	t['Fint{0}'.format(int(round(nu_freq)))] = [10**t['coefficients'][i][0] for i in range(len(t))]
	
	# Subsetting the new coefficients for testing purposes.
	t_co = np.array(t['coefficients'])
	
	#"""
	# Testing by plots:

	N = 10
	rand_ind = np.random.randint(0,len(t),N)

	nu_smooth = np.linspace(72.0,1400.0,1000)

	nu_smth_old = nu_smooth/old_freq
	nu_smth_nu = nu_smooth/nu_freq

	for j in rand_ind:

		Flux_old = 10**poly_calc(np.log10(nu_smth_old),t_co_old[j][::-1]).flatten()
		Flux_nu = 10**poly_calc(np.log10(nu_smth_nu),t_co[j][::-1]).flatten()

		print t_co_old[j]
		print t_co[j]

		plt.loglog(nu_smooth,Flux_old,label='Freq = {0} MHz'.format(int(round(old_freq))))
		#plt.semilogy(nu_smth_old,Flux_old,label='Freq = {0} MHz'.format(int(round(old_freq))))
		plt.loglog(nu_smooth,Flux_nu,label='Freq = {0} MHz'.format(int(round(nu_freq))),ls='--')
		#plt.semilogy(nu_smth_nu,Flux_nu,label='Freq = {0} MHz'.format(int(round(nu_freq))),ls='--')
		plt.xlabel(r'$\rm{\nu\,MHz}$',fontsize=18)
		plt.ylabel(r'$\rm{S\,Jy}$',fontsize=18)
		plt.legend()
		plt.show()
	#"""
	
	#t.write("Total-{0}MHz-skymodel.fits".format(int(round(nu_freq))),"w")
	
	if options.nu_catalogue != None:
		t.write("{0}.fits".format(str(options.nu_catalogue)),"w")
	else:
		t.write("Total-{0}MHz-skymodel.fits".format(str(options.nu_freq)),"w")

	end0 = time.time()

	print "Run time = {0} s".format(end0-start0)
