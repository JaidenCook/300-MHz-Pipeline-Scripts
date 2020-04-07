#!/usr/bin/python

__author__ = "Jaiden Cook"
__credits__ = ["Jaiden Cook"]
__version__ = "0.0.1"
__maintainer__ = "Jaiden Cook"
__email__ = "Jaiden.Cook@student.curtin.edu"

# Generic stuff:
import os,sys
import time
from datetime import datetime
import glob
import shutil
import re
from math import pi
import warnings
import subprocess
warnings.filterwarnings("ignore")

# Array stuff:
import numpy as np
warnings.simplefilter('ignore', np.RankWarning)

# Plotting stuff:
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
plt.style.use('seaborn-white')
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# Parser options:
from optparse import OptionParser

# Multiprocessing stuff:
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

# Scipy stuff:
import scipy
from scipy import stats

# Astropy stuff:
from astropy import wcs
from astropy.io import fits
from astropy.io import ascii
from astropy.io.votable import parse_single_table
from astropy.table import Table,Column,vstack
from astropy.io.votable import writeto as writetoVO

# This is a temporary hack.
sys.path.append(os.path.abspath("/home/jaiden/Documents/Masters_Project/bin"))
from JOOF import *

def subprocess_cmd(command,silent=False):
	process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	proc_stdout = process.communicate()[0].strip()
	if (silent==False):
		print proc_stdout
	return (proc_stdout)

if __name__ == "__main__":

	# Defining the parser:
	usage="Usage: %prog [options]\n"
	parser = OptionParser(usage=usage)
	parser.add_option('--freq_cent',dest="freq_cent",default=None,help="Input the normalisation frequency in MHz.\n")
	parser.add_option('--catalogue',dest="catalogue",default=None,help="Input list of sources with RA and Dec.\n")
	parser.add_option('--analysis',dest="analysis_cond",default=False,help="Output error and variance tables.\n")
	parser.add_option('--plot',dest="plot_cond",default=False,help="Save png plots of random SED's from the input catalogue..\n")
	parser.add_option('--filter_cat',dest="filter_cat",default=None,help="List of sources to filter from the main catalogue.\n")

	(options, args) = parser.parse_args()

	start0 = time.time()

	if options.catalogue.split('.')[1] == 'fits':
		# Case where the input file is a fits file.
		# Will cases for votables as well at a later time.
		# This code is for opening fits files:
		header = fits.getheader("{0}".format(options.catalogue))
		data = fits.getdata("{0}".format(options.catalogue))
	elif options.catalogue.split('.')[1] == 'csv':
		# Situation where the input catalogue is a csv file.
		# This is not recommended.
		# Useful for testing purposes
		data = ascii.read(options.catalogue)
	else:
		print "File format {0} not recognized exiting python!".format(options.catalogue.split('.')[1])
		sys.exit(0)

	
	# Table formats the fits file into a user readable output.
	t = Table(data) # Puts fits data into a readable table format.
	
	# For testing purposes:
	#t = t[:1000]
	
	col_names = t.colnames
	col_names_str = str(col_names)

	# Flux and error column names:
	flux_col_names = [string[string.find("S"):-1] for string in col_names_str[col_names_str.find("S"):col_names_str.find("e_S")].split(',')][:-1]
	e_flux_col_names = [string[string.find("e_S"):-1].strip("'") for string in col_names_str[col_names_str.find("e_S"):].split(',')]

	# Number of sources in the catalogue:
	N_sources = len(t[flux_col_names[0]])
	print "Number of sources in the catalogue = {0}".format(N_sources)

	# User inputed central frequency:
	freq_cent = float(options.freq_cent) # MHz
	print "Central frequency = {0} MHz".format(freq_cent)

	#### If filter catalogue is not None ####

	if options.filter_cat != None:
	
		if options.filter_cat.split('.')[1] == 'fits':
		
			# Case where the input file is a fits file.
			t_filt = Table(fits.getdata("{0}".format(options.filter_cat)))
		elif options.filter_cat.split('.')[1] == 'csv':
		
			# Situation where the input catalogue is a csv file.
			t_filt = ascii.read(options.filter_cat)
		else:
			"Filter catalogue format not recognized!"
			pass

		# RA and DEC filter source vectors:
		RA_filt_vec = t_filt["RA"]
		DEC_filt_vec = t_filt["DEC"]

		# Setting the angular separation threshold:
		ang_sep_thresh = 2.0/60.0 # in degrees.

		print "Filtering sources out of the main catalogue."

		for i in range(len(RA_filt_vec)):

			# Calculating the angular seperation between filter source and every other source:
			del_R_temp_vec = np.sqrt((np.array(t["RA"]) - RA_filt_vec[i])**2 + (np.array(t["DEC"]) - DEC_filt_vec[i])**2)

			#print "min(Delta R) = {0}".format(np.min(del_R_temp_vec))

			# Taking only the sources above the threshold:
			t_tmp = t[del_R_temp_vec < ang_sep_thresh]
			t = t[del_R_temp_vec > ang_sep_thresh]

			if i == 0:
				t_filt_puma = t_tmp
			else:
				t_filt_puma = vstack([t_filt_puma,t_tmp])

			# In future I may add an option to save the new thresholded catalouge.

		t_filt_puma.write("filtered_sources.fits","w")

		print "Total number of sources filtered = {0}".format(N_sources - len(t))
		# Setting the new number of sources after filtering.
		N_sources = len(t[flux_col_names[0]])
	else:
		pass

	# Initialising flux and error arrays, these should have the same dimensions:
	Flux_array = np.nan_to_num(np.array([np.array(t[chan]) for chan in flux_col_names]))
	e_Flux_array = np.nan_to_num(np.array([np.array(t[e_chan]) for e_chan in e_flux_col_names]))

	# Frequency list:
	Nu = np.array([float(flux_name.strip("S")) for flux_name in flux_col_names])

	################################################
	#"""
	# Note: this is not a good way to threshold problem sources.
	#
	if options.analysis_cond == "True":
		# Looking for odd sources in GLEAM:	
		GLEAM_flux_array = np.array(Flux_array[:20,:])
		GLEAM_flux_array[GLEAM_flux_array == 0.0] = np.nan
	
		GLEAM_chan_rms = scipy.stats.iqr(GLEAM_flux_array,axis=0,nan_policy='omit')

		# Might need to change this section.
		# Getting only the non-nan columns:
		trunc_GLEAM_flux_arr = GLEAM_flux_array[:,np.isnan(GLEAM_chan_rms) != True]
		t_var = t[np.isnan(GLEAM_chan_rms) != True]
		GLEAM_chan_rms = GLEAM_chan_rms[np.isnan(GLEAM_chan_rms) != True]
	
		# Setting the variable indx array:
		Var_ind_vec = (np.mean(trunc_GLEAM_flux_arr,axis=0) - 10*GLEAM_chan_rms) > 0.0
	
		# Subsetting the variable table:
		t_var = t_var[Var_ind_vec]
	
		print "Number of variable sources = {0}".format(len(t_var))
		print "Writing pumav5_var_table.fits table!"
		t_var.write("pumav5_var_table.fits","w")
	else:
		pass
	#"""
	################################################

	print "Frequency samples = {0} [MHz]".format(Nu)

	# Random sample for plotting purposes:
	#rand_ind_arr = np.arange(len(Flux_array[0,:]))

	# k=number of parameters per polynomial, which is order = k-1.
	# Specifying vector:
	#k_vec = np.array([2,3,4,5,6,7,8,9]) # Note: Don't exceed k=4, things get whacky.
	k_vec = np.array([2,3]) # Note: Don't exceed k=4, things get whacky.
	#k_vec = np.array([2]) # Note: Don't exceed k=4, things get whacky.

	print "Max order fitting = {0}".format(max(k_vec)-1)

	# Determining the polynomial coefficient array and the associated chisquared array:
	poly_co_array,Chisqd_arr = poly_coefficients(np.log10(Nu/freq_cent),np.log10(Flux_array),k_vec,(np.log(10)*Flux_array)/e_Flux_array,True) 
	
	# Specifying how many non-zero channels there are per source. concatenating into a vector:
	N_vec = np.array([len(Flux_array[:,i][Flux_array[:,i] > 0.0]) for i in range(N_sources)])

	# Using the BIC to determine the optimal order polynomial for each source.
	dBIC_thresh = 6.0
	Order_vec = BOOF(np.transpose(Chisqd_arr),k_vec,N_vec,dBIC_thresh)

	print "Delta BIC significance condition = {0}".format(dBIC_thresh)

	if options.analysis_cond == "True":

		# Determining the SI avg and sd:
		SI_avg = np.mean(poly_co_array[0,:,-2]) 
		SI_sd = scipy.stats.iqr(poly_co_array[0,:,-2])
		
		# Finding sources which are 5-sigma away from the mean:
		SI_index = np.logical_or(poly_co_array[0,:,-2] > SI_avg+5*SI_sd, poly_co_array[0,:,-2] < SI_avg-5*SI_sd)
	
		# Printing out the statistics:
		print "Number of sources SI > 5-sigma = {0}".format(len(poly_co_array[0,:,-2][SI_index]))
		print "Mean spectral index = {0}".format(np.round(SI_avg,4))
		print "Std spectral index ={0}".format(np.round(SI_sd,4))
		print "max(SI) = {0}".format(np.max(poly_co_array[0,:,-2]))
		print "min(SI) = {0}".format(np.min(poly_co_array[0,:,-2]))
	
		# Plotting the spectral index distribution:
		fig, axs = plt.subplots(1, figsize = (11,7.5), dpi=90, tight_layout=True)

		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		SI_txt_box = '\n'.join((r'$\mu=%.4f$' % (SI_avg, ),r'$\sigma=%.4f$' % (SI_sd, )))
		
		axs.hist(poly_co_array[0,:,-2], bins = 100)
		axs.set_xlabel(r'$\rm{SI}$', fontsize = 14)
		axs.set_ylabel(r'$\rm{N}$', fontsize = 14)
		axs.text(0.05, 0.95, SI_txt_box, transform=axs.transAxes, fontsize=14, verticalalignment='top', bbox=props)

		plt.savefig("SI-hist.png")

		# Diagnostic, looking for sources without positive spectral indices.
		t_err = t[SI_index]
	
		if len(t_err) > 1:
			print "Outputting table with positive spectral indices."
			t_err.write("pumav5_errtable.fits","w")
			print "pumav5_errtable.fits written!"
	else:
		pass

	print "Using the BIC to determine the optimal order fit to each source!"
	# There might be a quicker way to do this.
	# Might be worth parellizing this process if I can't speed it up with an array operation:
	Opt_ply_co_arr = np.array([poly_co_array[k_vec == (Order_vec+1)[i],i,:].flatten() for i in range(N_sources)])

	# Setting the maximum number of entries, this is for A-team sources, may increase this to 12.
	#Opt_ply_co_arr = np.insert(Opt_ply_co_arr,0,[[0.0],[0.0]],axis=1)
	Opt_ply_co_arr = np.insert(Opt_ply_co_arr,0,[[0.0],[0.0],[0.0]],axis=1)

	# Determining how many sources are first or second order fits:
	N_par = np.array([len(Opt_ply_co_arr[i,:][np.abs(Opt_ply_co_arr[i,:]) > 0.0]) for i in range(len(Opt_ply_co_arr[:,0]))])

	print "Number of first order fits = {0}".format(len(N_par[N_par == 2]))
	print "Number of second order fits = {0}".format(len(N_par[N_par == 3]))

	# Calculating the 300 MHz flux density.
	Fint = 10**poly_calc(0.0,Opt_ply_co_arr)

	print "Mean {2} MHz Flux density = {0} +- {1} [Jy]".format(np.round(np.mean(Fint),4),np.round(scipy.stats.iqr(Fint),4),freq_cent)
	print "Mean Spectral index = {0} +- {1}".format(np.round(np.mean(Opt_ply_co_arr[-1]),4),np.round(scipy.stats.iqr(Opt_ply_co_arr[-1])))

	#"""
	############This section of code is for testing purposes##############
	# Current working directory:

	if options.plot_cond == "True":
		
		# Setting the current working directory:
		pwd = subprocess_cmd("pwd",silent=False)
		
		# Setting the random index array:
		if len(Flux_array[0,:]) < 20:
			rand_ind_arr = np.random.randint(0,len(Flux_array[0,:]),len(Flux_array[0,:]))
		else:
			rand_ind_arr = np.random.randint(0,len(Flux_array[0,:]),200)
	
		# Make plot directory:
		subprocess_cmd("mkdir {0}/Fit_plots".format(pwd),silent=False)
	
		rand_ind_arr = np.array([300658,303692])
		# Generating plots:
		for i in rand_ind_arr:
			poly_fit_plot(Nu/freq_cent,Flux_array[:,i],poly_co_array[:,i,:],k_vec,e_Flux_array[:,i],plot=False)
			subprocess_cmd("mv Poly_fit.png {0}/Fit_plots/{1}Poly_fit.png".format(pwd,i),silent=False)
	else:
		pass

	Coeff_list = [tuple(Opt_ply_co_arr[i,::-1]) for i in range(N_sources)]


	print "Number of sources = {0}".format(len(t))


	print "Creating new table!"
	t_new = Table([np.array(t[col_names[0]]),np.array(t[col_names[1]]),np.array(t[col_names[2]]),np.array(t[col_names[3]]),np.array(t[col_names[4]])\
		,np.array(t[col_names[5]]), list(Fint), Coeff_list],names=col_names[:6]+["Fint{0}".format(int(freq_cent)),"coefficients"], meta={'name': 'first table'})
	
	print "Number of sources = {0}".format(len(t_new))

	print "Writing output table!"
	#t_new.write("pumav0_{0}Mhz_skymodel.fits".format(int(freq_cent)),"w")
	t_new.write("pumav1_{0}Mhz_skymodel.fits".format(int(freq_cent)),"w")

	if options.analysis_cond == "True":
		
		# Subsetting and writing a new bright sources catalogue:
		t_bright = t[np.argsort(Fint,axis=0)][len(t)-10:]
		t_bright.write("pumav5_bright.fits","w")

	else:
		pass

	print 'Completed'
	end0 = time.time()
	
	print "Total runtime = ", np.round(end0 - start0,3),"s\n"





# Filter sources before calling the catalogue into Sky_model.

"""
	#flux_col_names = [string[string.find("S"):-1] for string in col_names_str[col_names_str.find("S"):col_names_str.find("e_S")].split(',')][:-9]
	rel_err_list = np.array([np.mean(t[e_flux_col_names[i]][np.isnan(t[e_flux_col_names[i]]) == False]/t[flux_col_names[i]][np.isnan(t[e_flux_col_names[i]]) == False]) for i in range(len(flux_col_names[:-8]))])
	iqr_rel_err_list = np.array([scipy.stats.iqr(t[e_flux_col_names[i]][np.isnan(t[e_flux_col_names[i]]) == False]/t[flux_col_names[i]][np.isnan(t[e_flux_col_names[i]]) == False]) for i in range(len(flux_col_names[:-8]))])



	# Janky, jank jank jank
	#flux_col_names = [string[string.find("S"):-1] for string in col_names_str[col_names_str.find("S"):col_names_str.find("e_S")].split(',')][-1]
	#e_flux_col_names = [string[string.find("e_S"):-1].strip("'") for string in col_names_str[col_names_str.find("e_S"):].split(',')]




	cent_chan = float(int(np.mean(Freq_list)))
	print "Channels = {0}".format(Freq_list)
	print "Middle channel = {0}".format(cent_chan)

	X = np.linspace(min(Freq_list),max(Freq_list),1000)/cent_chan
	
	fig, ax = plt.subplots(1)
	ax.errorbar(Freq_list/cent_chan,rel_err_list,yerr=iqr_rel_err_list,ecolor='k',fmt='o',markeredgecolor='k',capsize=1.5,elinewidth=0.5,label='Relative Error')
	#ax.scatter(Freq_list/cent_chan,rel_err_list,edgecolor='k')
	ax.fill_between(Freq_list/cent_chan, rel_err_list - iqr_rel_err_list, rel_err_list + iqr_rel_err_list, facecolor='blue', alpha=0.2,label=r'$\rm{1-\sigma}$')

	# JOOF.py might be useful here.
	poly_co_array = poly_coefficients(Freq_list/cent_chan,rel_err_list,k_max=len(Freq_list)-1,weights=None)#,log=False)

	#ax.fill_between(X, np.polyval(poly_co_array[2,:],X) - np.mean(iqr_rel_err_list), np.polyval(poly_co_array[2,:],X) + np.mean(iqr_rel_err_list), facecolor='blue', alpha=0.2,label='1 sigma range')
	ax.legend(loc='upper left')


	Chi_squared = np.array([np.sum(((rel_err_list - np.polyval(poly_co_array[i,:],Freq_list/cent_chan))/iqr_rel_err_list)**2) for i in range(len(poly_co_array[0,:]))])
	BIC = np.array([Chi_squared[i] + len(poly_co_array[i,:][poly_co_array[i,:] > 0])*np.log(len(Freq_list)) for i in range(len(poly_co_array[0,:]))])

	min_order = np.argmin(BIC)

	print "Minimum order = {0}".format(min_order)

	# Defining lambda function this will be useful later.
	poly_err_cal = lambda p,min_ord,x0,x: np.polyval(p[min_ord,:],x/x0)

	print poly_err_cal(poly_co_array,min_order,cent_chan,170.0)
	print poly_err_cal(poly_co_array,min_order,cent_chan,190.0)
	print poly_err_cal(poly_co_array,min_order,cent_chan,210.0)


	ax.errorbar(170.0/cent_chan,np.polyval(poly_co_array[min_order,:],170.0/cent_chan),markeredgecolor='k',fmt='o',color='red',label='S170')
	ax.errorbar(190.0/cent_chan,np.polyval(poly_co_array[min_order,:],190.0/cent_chan),markeredgecolor='k',fmt='o',color='yellow',label='S190')
	ax.errorbar(210.0/cent_chan,np.polyval(poly_co_array[min_order,:],210.0/cent_chan),markeredgecolor='k',fmt='o',color='purple',label='S210')

	for i in range(len(poly_co_array[0,:])):

		if i % 2 == 0 and i < 10 and i > 0:
			plt.plot(X,np.polyval(poly_co_array[i,:],X),label="Polynomial Order = {0}".format(i))

	ax.set_ylabel(r'$\rm{\langle \frac{u(S)}{S} \rangle}$',fontsize = 14)
	ax.set_xlabel(r'$\rm{\nu/\nu_o\:[MHz]}$',fontsize = 14)
	plt.tight_layout()
	plt.legend()
	#plt.show()


	t['e_S170'] = 0.0
	t['e_S190'] = 0.0
	t['e_S210'] = 0.0

	t['e_S170'][t['S170'] > 0 ] = t['S170'][t['S170'] > 0 ]*np.polyval(poly_co_array[min_order,:],170.0/cent_chan)
	t['e_S190'][t['S190'] > 0 ] = t['S190'][t['S190'] > 0 ]*np.polyval(poly_co_array[min_order,:],190.0/cent_chan)
	t['e_S210'][t['S210'] > 0 ] = t['S210'][t['S210'] > 0 ]*np.polyval(poly_co_array[min_order,:],210.0/cent_chan)

	# After adding the error columns we have to redefine col_names, and e_flux_col_names.
	col_names = t.colnames
	col_names_str = str(col_names)

	e_flux_col_names = [string[string.find("e_S"):-1].strip("'") for string in col_names_str[col_names_str.find("e_S"):].split(',')]

	print e_flux_col_names

	# Initialising flux and error arrays:
	Flux_array = np.nan_to_num(np.array([np.array(t[chan]) for chan in flux_col_names]))
	e_Flux_array = np.nan_to_num(np.array([np.array(t[e_chan]) for e_chan in e_flux_col_names]))
	
	print "shape(Flux_array) = {0}".format(np.shape(Flux_array))
	print "shape(e_Flux_array) = {0}".format(np.shape(e_Flux_array))

	# Get rid of sources which only have one flux density measurement, polyfit takes a vector as an input:
	ind_array = np.array([len(Flux_array[:,i][Flux_array[:,i] > 0.0]) for i in range(N_sources)])#####
		
	print "Flux array shape before filtering; shape(Flux_array) = {0}".format(np.shape(Flux_array))
	Flux_array = Flux_array[:,ind_array > 1]
	e_Flux_array = e_Flux_array[:,ind_array > 1]
	print "Flux array shape after filtering; shape(Flux_array) = {0}".format(np.shape(Flux_array))

	# Set flux density values that are negative and their errors to be zero.
	e_Flux_array[Flux_array < 0] = 0.0 # Error array has to come first.
	e_Flux_array[e_Flux_array < 0] = 0.0 # Setting negative errors as zero.
	Flux_array[Flux_array < 0] = 0.0 # Setting negative fluxes as zero.

	# Specifying an array of the relative errors that has the same dimensions as the table.
	# This way I can do a simple matrix multiplication.
	rel_err_tot_arr = np.ones([len(flux_col_names),len(Flux_array[0,:])])
	rel_err_tot_arr[:len(rel_err_list),:] = np.transpose(rel_err_list*np.transpose(rel_err_tot_arr[:len(rel_err_list),:]))
	rel_err_tot_arr[rel_err_tot_arr == 1] = 0.1
	
	# Specifying an index array where there are fluxes without errors:
	index_array = np.logical_and(Flux_array > 0.0,np.abs(e_Flux_array) == 0.0)#####

	# Setting the estimated errors:
	e_Flux_array[index_array] = Flux_array[index_array]*rel_err_tot_arr[index_array]
	
	# Determining the number of sources that have missing errors and fluxes.
	N_flux = np.array([len(Flux_array[:,i][Flux_array[:,i] > 0.0]) for i in range(len(Flux_array[0,:]))])
	N_e_flux = np.array([len(e_Flux_array[:,i][e_Flux_array[:,i] > 0.0]) for i in range(len(e_Flux_array[0,:]))])

	print "Number of sources with more flux channels that errors = {0}".format(len(N_flux[(N_flux - N_e_flux) > 0.0]))
	print "Number of sources with more error flux channels that errors = {0}".format(len(N_flux[(N_flux - N_e_flux) < 0.0]))
	print "Max number of missing flux channels = {0}, index = {1}".format(np.max(N_flux - N_e_flux),np.argmax(N_flux - N_e_flux))	
	print "Max number of missing error flux channels = {0}, index = {1}".format(-1*np.min(N_flux - N_e_flux), np.argmin(N_flux - N_e_flux))
	
	
	# Need to index these columns.
	Other_cols = list(zip(*list(Table.as_array(t[col_names[:7]][ind_array > 1])))) 

	flux_list = [list(Flux_array[i,:]) for i in range(len(Flux_array[:,0]))]
	e_flux_list = [list(e_Flux_array[i,:]) for i in range(len(e_Flux_array[:,0]))]

	proto_table = Other_cols + flux_list + e_flux_list

	print np.shape(proto_table)
	print len(col_names)

	t_new = Table(proto_table, names=col_names, meta={'name': 'first table'})

	print t_new

	t_new.write("pumav5_full_spectral.fits","w")
"""