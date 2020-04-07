#!/usr/bin/python

__author__ = "Jaiden Cook"
__credits__ = ["Jaiden Cook"]
__version__ = "1.1.2"
__maintainer__ = "Jaiden Cook"
__email__ = "Jaiden.Cook@student.curtin.edu"

# Standard library:
import os
import sys
import time
from datetime import datetime
import glob
import shutil
import re
from math import pi
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
plt.style.use('seaborn-white')
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# Array stuff:
import numpy as np
warnings.simplefilter('ignore', np.RankWarning)

# Parser options:
from optparse import OptionParser

# Multiprocessing stuff:
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

# Scipy stuff:
from scipy.special import binom
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


def poly_calc(x,coeff_arr):
	"""
	Takes an input coefficient array and calculates the polynomial for N data points, for M sources.
	Outputs a 2D array, where each row is a source, each column is the polynomial value(s) at the input
	x value(s).

	Args:
	x : float/vector; Arbitrary x variable, when log=True this should be normalised by some constant x-value
	to ensure the argument of the logarithm is dimensionless.
	coeff_arr : array_like; This is an array or vector of polynomial coefficients
	"""

	# Converting single value variables to an array of length 1. For use in Numpy.
	x = x*np.ones(1)

	if len(np.shape(coeff_arr)) > 1:
		
		poly_arr = np.zeros([len(coeff_arr[:,0]),len(x)])
		k = len(coeff_arr[0,:])

	else:
		# For cases when you only want to calculate for a single source.
		poly_arr = np.zeros(len(x))
		k = len(coeff_arr)

	order = k - 1
	
	for i in range(k):

		if len(np.shape(coeff_arr)) > 1:
			
			poly_arr = poly_arr + np.outer(coeff_arr[:,i],(x)**(order - i))
		
		else:
			# For cases when you only want to calculate for a single point.
			poly_arr = poly_arr + np.outer(coeff_arr[i],(x)**(order - i))

	return poly_arr

def poly_coefficients(x,y,k_vec,weights=None,full_cond=False):
	"""
	This function performs a parallel np.polyfit() on a 2D array, where the number of rows
	is the number of sources to fit. The order is generalised, and the user has the option
	of using a log or normal polynomial. This function returns a 3D array where each initial 
	entry corresponds to the order, the 2D array for that order contains the polynomial 
	coefficients for each row (source). For orders less than the maximum order, additional
	column entries are set to 0.

	Args:
	x : vector; This is a vector of x-positions/values which correspond to the input y-values.
	y : array_like; This is a 2D array of y-values where each row is a different source. Future versions will
	allow for single vectors of y-values.
	k_max : vector; This is the maximum number of polynomial parameters polyfit will fit to. The order
	of the polynomial is defined as order = k - 1 for an arbitrary polynomial, hence 
	max{order} = k_max - 1.
	weights : vector; Default is None, for Gaussian weights this takes the form of 1/Sigma.
	log = False : bool; This function can work with log-polynomials and normal space polynomials, the
	default is set to normal space.
	"""

	# Maximum number of parameters:
	k_max = np.max(np.array(k_vec)) # Wrapped in an array.

	# Be careful with the dimensions of the input y array.
	if len(np.shape(y)) == 1:

		# Case in which we have a 1D input array.
		poly_co_arr = np.zeros([len(np.array(k_vec)),k_max])
	
		if cov_cond != False:
			Chisqd_vec = np.zeros(len(np.array(k_vec)))

		i = 0 # Counting variable for indexing purposes.
		for k in k_vec:

			print "Fitting polynomial order = {0}".format(k-1)

			poly_co_vec = np.polyfit(x,y,k,w=weights,full=full_cond)

			if full_cond != False:
				# Case when the you want the chisqd value
				poly_co_arr[i,(k_max-1)-k+1:] = poly_co_vec[0]
				Chisqd_vec[i] = poly_co_vec[1]
			else:
				poly_co_arr[i,(k_max-1)-k+1:] = poly_co_vec
			
			i += 1

		if full_cond != False:
			
			return poly_co_arr, Chisqd_vec
		
		elif full_cond == False:
			
			return poly_co_arr
		else:
			return None

	elif len(np.shape(y)) == 2:

		# Case in which we have a 2D input array.
		N_sources = len(y[0,:])
		poly_co_arr = np.zeros([len(np.array(k_vec)),N_sources,k_max])
		
		if full_cond != False:
			Chisqd_array = np.zeros([len(np.array(k_vec)),N_sources])

		i = 0 # Counting variable for indexing purposes.
		for k in k_vec:
			
			print "Fitting polynomial order = {0}".format(k-1)
			
			if np.any(weights) != None:
			
				# Case when you want to fit with weights.
				# is nan case when dealing with polynomials which have zero values.
				# is inf case deals with log-polynomials, where log(0) = -inf.
				poly_co_vec = np.array(Parallel(n_jobs=-1)(delayed(np.polyfit)(x[np.isnan(y[:,j])!=True][np.isinf(y[:,j])!=True],\
					y[:,j][np.isnan(y[:,j])!=True][np.isinf(y[:,j])!=True],k-1,w=weights[:,j][np.isnan(y[:,j])!=True][np.isinf(y[:,j])!=True],full=full_cond)\
					 for j in tqdm(range(N_sources))))
			else:
				
				# Case when you want an unweighted fit.
				poly_co_vec = np.array(Parallel(n_jobs=-1)(delayed(np.polyfit)(x[np.isnan(y[:,j])!=True][np.isinf(y[:,j])!=True],\
					y[:,j][np.isnan(y[:,j])!=True][np.isinf(y[:,j])!=True],k-1,w=weights,full=full_cond) for j in tqdm(range(N_sources))))

			if full_cond != False:
				# Case when the you want the chisqd value
				poly_co_arr[i,:,(k_max-1)-k+1:] = np.array(list(poly_co_vec[:,0]))

				# Index source which have Chi-squared values, a value of 0 means that N = k for that source.
				temp_ind_arr = np.array([len(y[:,l][np.isnan(y[:,l])!=True][np.isinf(y[:,l])!=True]) for l in range(N_sources)]) > k
				
				# Inputing chis squared values into the chisquared array:
				Chisqd_array[i,:][temp_ind_arr] = poly_co_vec[temp_ind_arr,1]
			else:
				poly_co_arr[i,:,(k_max-1)-k+1:] = poly_co_vec
			
			i += 1

		if full_cond != False:
			
			return poly_co_arr, Chisqd_array
		
		elif full_cond == False:
			
			return poly_co_arr
		else:
			return None

def log_poly_co_map(a,x_a,x_b):
	"""
	This function takes the log-polynomial coefficients for some nth order polynomial, normalised at
	some value x_a, and maps them to the coefficients of the same log-polynomial normalised at some
	value x_b. This has been generalised to work for one polynomial fit to one set of data, or for
	N_sources polynomials fitted to N_sources of data, in this case the input is a 2D array, where
	each row is a set of coefficients fitted to some set of data.
 
	Args:
	a : vector_like, array_like; Polynomial coefficients for order n polynomial a = [a0,a1,a2,...,an],
	or a = [[a0,0 a0,1 ... a0,n],...,[an,0 an,1 ... an,n]].
	x_a : float_like; Normalisation value for polynomial p_a.
	x_b : float_like; Normalisation value for polynomial p_b.
	"""

	# The for loops can be eliminated and replaced with arrays, and outer product operations.

	if len(np.shape(a)) == 1:

		# Vector case:
		n = len(a) # Number of parameters.
		b = np.zeros(n) # Initialising p_b coefficient list.
	
		# Generating function:
		for j in range(n):
			element = 0
			for i in range(j,n):
				element = element + binom(i,i-j)*a[i]*((np.log10(x_b/x_a))**(i-j))
	
			b[j] = element
	
		return b

	elif len(np.shape(a)) == 2:

		# 2D array case:
		n = len(a[0,:]) # Number of parameters.
		N_sources = len(a[:,0])
		b = np.zeros(np.shape(a)) # Initialising p_b coefficient list.
	
		# Generating function:
		for j in range(n):
			element = np.zeros(N_sources)
			for i in range(j,n):
				element = element + binom(i,i-j)*a[:,i]*((np.log10(x_b/x_a))**(i-j))
	
			b[j,:] = element
	
		return b
	
	else:
		# In the event a 3D or higher dimensional array is given, return None:
		return None		


def poly_fit_plot(Nu_norm,Beam_vec,poly_co_list,k_vec,y_errors=None,Zen=None,theta=None,radec=False,plot=False):
	"""
	This function plots the log-beam response across the bandwidth, and the
	even or odd (depending on the user input) log-polynomial fits to the 
	log-beam response. The secondary plot is the squared residuals for each fitted
	order.

	Args:

	Nu_norm: Vector; Normalised frequecny vector, principally can be any generic x vector.
	Beam_vec: Vector; Beam value, principally can be any generic y vector.
	Poly_co_list: Array; Array of polynomial coefficients. Each row is a sucessive order, and each column is 
	the coefficient for that row in decending order, from the highest to lowerst coefficeint.
	y_errors: Vector; Vector of associated y-errors.
	k_vec: Vector; Vector of synonomous to order, whe
	Zen/Dec: Scalar, default = None; The Azimuth angle. Not None in the non-generic polynomial case. This is
	the declination when radec = True.	
	Theta/RA: Scalar, default = None; The Azimuth angle. Not None in the non-generic polynomial case. This is
	the RA when radec = True.
	even: Boolean, default = True; Plot only the even polynomials when True,plot odd polynomials when False.
	radec: Boolean, default = False; Plot local polar coords when False, plot RA/DEC coords when True. 
	plot: Boolean, default = False; If True, show plot, if False don't show plot.
	"""

	# Initialising:
	
	fig = plt.figure(figsize = (11,9), dpi=70, tight_layout=True) 
	gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
	
	# Adding top axis:
	ax0 = plt.subplot(gs[0], yscale='log',xscale='log')
	
	ax1 = fig.add_subplot(gs[1], sharex=ax0, yscale='log',xscale='log')
	
	#ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

	#ax0.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
	
	#plt.setp(ax0.get_xticklabels(), visible=False)	
	gs.update(hspace=0)

	# Necessary parameters:
	k_max = len(poly_co_list[:,0])

	# Determine how many parameters there are per entry:
	poly_k_vec = np.array([len(poly_co_list[i,:][np.abs(poly_co_list[i,:]) > 0.0]) for i in range(len(poly_co_list[:,0]))])

	# Creating an array of indices to plot:
	order_index_list = np.array([np.arange(len(poly_k_vec))[poly_k_vec == k_vec[i]] for i in range(len(k_vec))]).flatten()
	
	# Indexing for non-zero beam values:
	index_list = Beam_vec > 0.0

	# Plotting the initial log-beam response:
	if np.any(y_errors) != None:
		# Case where there are input y-errors.

		# Setting errors greater than y to be 0.999*y as to not break the plotting.
		y_errors_assym = y_errors
		y_errors_assym[(y_errors-Beam_vec) > 0.0] = 0.999*Beam_vec[(y_errors-Beam_vec) > 0.0]

		ax0.errorbar(300*Nu_norm[index_list],Beam_vec[index_list],yerr=[y_errors_assym[index_list],y_errors[index_list]]\
			,ecolor='k',fmt='o',markeredgecolor='k',capsize=1.5,elinewidth=0.5,label='Channels {0}-{1} [MHz]'.format(int(300*np.min(Nu_norm)),int(300*np.max(Nu_norm))))
		# Setting the y-limits
		#ax0.set_ylim([np.min(Beam_vec[index_list])/2.5,2.5*np.max(Beam_vec[index_list])])

		print np.min(Beam_vec[index_list])/2.5,2.5*np.max(Beam_vec[index_list])
	else:
		# Case where no errors are input.
		ax0.scatter(300*Nu_norm[index_list],Beam_vec[index_list],label=r'$\rm{Beam/Chan}$',edgecolors = 'k')
		# Setting the y-limits
		#ax0.set_ylim([np.min(Beam_vec[index_list])/1.1,1.1*np.max(Beam_vec[index_list])])
		#ax0.set_ylim([10e-5,2*10e-3])

	if radec == True and Zen != None and theta != None:
		# In this case, Zenith and Azimuth are replaced by Dec and Ra:
		plt.title("RA = {0} deg, DEC = {1} deg".format(np.round(np.degrees(theta),2)\
			,np.round(np.degrees(Zen),2)),fontsize = 14)
	elif radec == False and Zen != None and theta != None:
		# In this case were are dealing with local polar coordinates:
		plt.title("Azimuth = {0} deg, Zenith = {1} deg".format(np.round(np.degrees(theta),2)\
			,np.round(np.degrees(Zen),2)),fontsize = 14)

	# Creating an array of smooth points:
	Nu_norm_smth = np.linspace(np.min(Nu_norm[index_list]),np.max(Nu_norm[index_list]),1000)
	#ax1.set_xticks([0.26,0.5,4],minor=False)
	
	i = 0
	for index in order_index_list:

		# Setting the figure:
	
		# Fit axis:
		## Plotting the polynomial fits:
		ax0.plot(300*Nu_norm_smth,10**poly_calc(np.log10(Nu_norm_smth),poly_co_list[index,:])[0],ls='--',label="Polylog Order={0}".format(k_vec[i]-1))
		
		# Residual axis:
		ax1.scatter(300*Nu_norm[index_list],((10**poly_calc(np.log10(Nu_norm)[index_list],poly_co_list[index,:])[0] - Beam_vec[index_list])/Beam_vec[index_list])**2,\
			edgecolors = 'k',label='order={0}'.format(k_vec[i]-1))

		
		i += 1 # Counter for indexing the k_vec inside the for loop.

	ax0.set_ylabel(r'S [Jy]',fontsize = 20)
	ax1.set_ylabel(r'$\rm{\chi^2}$',fontsize = 20)
	

	# Setting plotting parameters:
	ax0.set_yscale('log')
	#ax0.set_yticks([np.round(np.min(10**poly_calc(np.log10(Nu_norm_smth),poly_co_list[index,:])[0]),2), \
	#	np.median(np.round((10**poly_calc(np.log10(Nu_norm_smth),poly_co_list[index,:])[0]),2)), \
	#	np.round(np.max(10**poly_calc(np.log10(Nu_norm_smth),poly_co_list[index,:])[0]),2)])#,minor=False)
	#ax0.get_yaxis().set_major_formatter(mtick.ScalarFormatter())


	ax1.set_xscale('log')
	#ax1.set_xticks([300*np.round(np.min(Nu_norm_smth),2)+20, 300*1.00, 300*np.round(np.max(Nu_norm_smth),2)-400],minor=False)
	ax1.set_xticks([300*np.round(np.min(Nu_norm_smth),2), 300*1.00, 300*np.round(np.max(Nu_norm_smth),2)])#,minor=False)
	ax1.get_xaxis().set_major_formatter(mtick.ScalarFormatter())
	
	
	ax0.yaxis.set_minor_locator(mtick.AutoMinorLocator(0.1))
	#ax1.yaxis.set_minor_locator(mtick.AutoMinorLocator(10))
	ax1.xaxis.set_minor_locator(mtick.AutoMinorLocator(10))

	ax0.tick_params(axis="y", labelsize=20,  width=1)
	ax1.tick_params(axis="x", labelsize=20,  width=1)
	
	ax0.legend(fontsize=20-i)
	
	ax1.set_xlabel(r'$\rm{\nu\,[MHz]}$', fontsize = 20)
	ax1.set_ylabel(r'$\rm{\chi^2}$', fontsize = 20)
	ax1.tick_params(axis="x", labelsize=20)
	ax1.tick_params(axis="y", labelsize=20)

	#plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=True)

	#plt.setp(ax0.get_yticklabels(), visible=True)
	#plt.setp(ax1.get_xticklabels(), visible=True)
	#plt.setp(ax1.get_yticklabels(), visible=True)
	# The figure is automatically saved:
	plt.savefig("Poly_fit.png")
	print "Poly_fit.png saved!"

	if plot == True:

		# Figure is only shown if specified:
		plt.show()
	else:
		plt.close()

	return None

def OOF_plot(k_max,OOF_ar,Zen=None,theta=None,plot=False):
	"""
	This function plots the Optimal Order Function (OOF) for a given source, which is identified
	by the index variable ind.

	Args:

	k_max: Scalar; the maximum number of free parameters.
	OOF_ar: Vector; The Optimal Order Function vector
	Zen: Scalar, default = None; The zenith angle. Not None in the non-generic polynomial case.
	Azimuth: Scalar, default = None; The Azimuth angle. Not None in the non-generic polynomial case.
	"""

	if Zen != None and theta != None:

		# 1D case
		plt.semilogy()
		plt.plot(range(k_max),OOF_ar,label ='OOF')
		plt.scatter(range(k_max),OOF_ar,edgecolor='k',label='Poly Min Order = {0}'.format(np.argmin(OOF_ar)))
		plt.xlabel('Order = k + 1',fontsize=14)
		plt.ylabel(r'$\rm{OOF(k)}$',fontsize=14)
		plt.title("Zenith = {0} deg, Azimuth = {1} deg".format(np.round(np.degrees(Zen),2),np.round(np.degrees(theta),2)),fontsize = 14)
		plt.legend()
		plt.tight_layout()
		plt.savefig("OOF_plt.png")
		print "OOF_plt.png saved!"
	
	elif Zen == None and theta == None:

		# 1D case no zenith or azimuth generic polynomial.
		plt.semilogy()
		plt.plot(range(k_max),OOF_ar,label ='OOF')
		plt.scatter(range(k_max),OOF_ar,edgecolor='k',label='Poly Min Order = {0}'.format(np.argmin(OOF_ar)))
		plt.xlabel('Order = k + 1',fontsize=18)
		plt.ylabel(r'$\rm{OOF(k)}$',fontsize=18)
		plt.legend(fontsize=18)
		plt.tight_layout()
		plt.savefig("OOF_plt.png")
		print "OOF_plt.png saved!"

	if plot == True:
		plt.show()
	else:
		plt.close()

	return None

def OOF_analysis_plots(Zen_samp,Az_samp,Ord_list,resid_list,plot=False):
	"""
	This function takes the output order and residuals from OOF, and generates a polar
	plot and histogram for both.

	Args:

	Zen_samp: Vector like; Sample of zenith values to be plotted.
	Az_samp: Vector like; Sample of azimuthal values to be plotted.
	Ord_list: Vector like; Corresponding optimal order values for each set of Zen/Az.
	resid_list: Vector like; Corresponding residual square summed values for each set of Zen/Az.
	plot: Boolean, default = False; If False dont plot, if true plot.

	Notes: Future versions of this will allow for grid plotting.
	"""	
	Ord_list = np.array(Ord_list)
	
	# Setting the text box parameters:
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	textstr_X = '\n'.join((r'mode=%.2f' % (np.float(stats.mode(Ord_list)[0]), ),r'$\mu=%.2f$' % (np.mean(Ord_list), )))
	textstr_Y = '\n'.join((r'$\mu=%.3e$' % (np.mean(resid_list), ),r'$\sigma=%.3e$' % (np.std(resid_list), )))

	r = np.sin(Zen_samp)
	
	fig = plt.figure(figsize = (12.5,9.5), dpi=90)
	
	# Plotting the output:
	ax1 = fig.add_subplot(221,projection="polar")
	ax1.set_title('Polar Order Map')
	c = ax1.scatter(Az_samp,r,c=Ord_list,cmap='viridis',label='Order')
	ax1.set_rmin(0.0)
	ax1.set_rmax(1.01)
	ax1.set_theta_zero_location('N')
	fig.colorbar(c,cmap='viridis')
	
	ax2 = fig.add_subplot(222, projection='polar')
	ax2.set_title('Polar Residual Map')
	c = ax2.scatter(Az_samp,r,c=np.log10(resid_list),cmap='viridis',label='Residuals')
	ax2.set_rmin(0.0)
	ax2.set_rmax(1.01)
	ax2.set_theta_zero_location('N')
	fig.colorbar(c,cmap='viridis')
	
	ax3 = fig.add_subplot(223)
	ax3.hist(Ord_list,bins=range(24),edgecolor='k',color='g',alpha=0.5,label='OOF k min',align='left')
	ax3.set_xlabel("Order",fontsize=18)
	ax3.set_title('Order Hist')
	ax3.text(0.75, 0.95, textstr_X, transform=ax3.transAxes, fontsize=14,verticalalignment='top', bbox=props)
	
	ax4 = fig.add_subplot(224)
	ax4.hist(resid_list,align='left',bins=25,range=(np.min(resid_list),2*np.mean(resid_list)),alpha=0.5,edgecolor='k',color='g',label='Residuals')
	ax4.set_xlabel("Residuals",fontsize=18)
	ax4.set_title('Residual Hist')
	ax4.text(0.65, 0.95, textstr_Y, transform=ax4.transAxes, fontsize=14,verticalalignment='top', bbox=props)
	ax4.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
	
	plt.tight_layout()
	plt.savefig("OOF_plot.png")
	print "OOF_plot.png saved!"
	
	if plot == True:
		plt.show()
	else:
		plt.close()

	return None

def OOF(Nu,Beam_arr,px_arr,figcond=False,Zen=None,theta=None,poly_co_list=None):
	"""
	The Optimal Order Function (OOF):
	This calculates the OOF value for each polynomial fit N_chans, by comparing the Residual
	Square Sum (RSS), and the number of squared channels N_chans, normalised by the squared
	degrees of freedom (N_chans - k = dof). These two parameters are normalised, and the 
	OOF is determined by taking the magnitude of a vector described in the parameter space of
	these two functions.

	Args:
	Nu : array_like; Vector of frequency values in Hz, this should be normalised by the central frequency 
	of the fitting bandwidth.
	Beam_arr : array_like; 2D array of MWA coarse channel beam values for N_sources.
	px_arr : array_like; 3D array of polynomial MWA coarse channel beam approximations across the band with
	for k_max - 1 orders, for each source.
	figcond = False : bool; Save analysis figures; random sample of the log-beam response fit with 
	k_max - 1 order polynomials across the bandwidth for a single zenith  and azimuth angle;
	the OOF of the same random sample; and the overal statistical histograms and polar plots.
	Zen = None : vector; The zenith angle in radians for each pointing in the sky.
	theta = None : vector; The azimuth angle in radians for each pointing in the sky.
	poly_co_list = None : array_like; List/array of polynomial coefficients fit for each source, necessary 
	for calling nested plotting functions.
	"""

	if len(np.shape(px_arr)) == 2:

		# This is the case where we only want to operate on a single source.

		# Specifying necessary parameters:
		k_max = len(px_arr[:,0])# Should be fixed to the array dimensions.
		N_chans = np.float(len(Nu))
		k_vec = (np.arange(k_max) + 1).astype(float) # Number of parameters
		dof = N_chans - k_vec
		plotcond = False # Plot figures yes = True, no = False.

		## Determining the residuals for each channel, for each polynomial fit.
		RS = ((px_arr - Beam_arr)/Beam_arr)**2 
		RSS = np.sum(RS,axis=1)

		# Number of points normalised by the degrees of freedom:
		N_norm_dof = (N_chans/dof)**2
		N_norm_max = (N_chans)**2 # This is the max function value.
	
		# Calculating OOF:
		## Normalised to first order, since we don't care about the zeroth order.
		OOF = np.sqrt((N_norm_dof/N_norm_max)**2 + (RSS/RSS[1])**2)/np.sqrt((N_norm_dof[1]/N_norm_max)**2  + 1)
	
		# Finding the minimum order, this corresponds to the index:
		Order_array = np.argmin(OOF, axis=0)

		if figcond == True:

			OOF_plot(k_max,OOF,plot=plotcond)

		elif figcond == False:
	
			# When cond = False return only the order array.
			return Order_array

	elif len(np.shape(px_arr)) > 2:

		# This is the case where we want to operate on an array of sources.

		# Specifying necessary parameters:
		k_max = len(px_arr[:,0,0])# Should be fixed to the array dimensions.
		N_sources = len(px_arr[0,:,0])
		N_chans = np.float(len(Nu))
		k_vec = (np.arange(k_max) + 1).astype(float) # Number of parameters
		dof = N_chans - k_vec
		plotcond = False # Plot figures yes = True, no = False.
	
		# Calculating the RSS:
		## dim(px_arr) = (k_max,N_sources,N_chans)
		## dim(Beam) = (N_sources,N_chans)
		##
		## Determining the residuals for each channel, for each polynomial fit.
		RS = ((px_arr - Beam_arr)/Beam_arr)**2 
		RSS = np.sum(RS,axis=2) # Summing the residuals for each channel.	
		##
		## dim(RSS) = (k_max,N_sources)
		##
	
		# Number of points normalised by the degrees of freedom:
		N_norm_dof = np.transpose(((N_chans/dof)**2)*np.ones([N_sources,k_max]))
		N_norm_max = (N_chans)**2 # This is the max function value.
	
		# Calculating OOF:
		## Normalised to first order, since we don't care about the zeroth order.
		OOF = np.sqrt((N_norm_dof/N_norm_max)**2 + (RSS/RSS[1,:])**2)/np.sqrt((N_norm_dof[1]/N_norm_max)**2  + 1)
	
		# Finding the minimum order, this corresponds to the index:
		Order_array = np.argmin(OOF, axis=0)
	
		if figcond == True:
	
			# Select this option when you want to use OOF_plot.
			Residual_array = np.array([RSS[:,i][Order_array[i]] for i in range(len(Order_array))])
	
			rand_ind = np.random.randint(0,N_sources)
	
			if np.any(Zen) != None and np.any(theta) != None:
	
				if np.any(poly_co_list[0]) != None:
					
					# Plotting the fitted log-polynomials to the log beam response.
					#Fit_plot(Nu_norm,Beam_vec,poly_co_list,y_errors=None,Zen=None,theta=None,even=True,radec=False,plot=False):
					poly_fit_plot(Nu,Beam_arr[rand_ind,:],poly_co_list[:,rand_ind,:],np.arange(2,14,2),None,Zen[rand_ind],theta[rand_ind],plot=plotcond)
				else:
					pass
	
				print "Plotting OOF analysis plots!"
				OOF_plot(k_max,OOF[:,rand_ind],Zen[rand_ind],theta[rand_ind],plot=plotcond)
				OOF_analysis_plots(Zen,theta,Order_array,Residual_array,plot=plotcond)
	
			# If no Zenith or Azimuth angle given, residual array is returned.
			return Order_array, Residual_array, OOF
	
		elif figcond == False:
	
			# When cond = False return only the order array.
			return Order_array

def BOOF(Chisqd_arr,k_vec,N_vec,dBIC_thresh,BICcond=False):
	"""
	The Baysian Information Criterion (BIC) Optimal Order Function (BOOF):
	This calculates the optimal order polynomial fit for N_sources, where N_sources is the length of the
	N_vec vector. The optial order is determined by calculating the BIC for each polynomial fit, which is
	given by:

	BIC = Chi2 + k*ln(n)

	where Chi2 is the chisquared or the residual squared sum, k is the number of parameters in the fit,
	and n is the number of data points that are fitted.

	Args:
	Chidqd_arr : array_like; Array which contains the chisquares value for each polynomial fit.
	k_vec : vector_like; Vector which contains the number of parameters per polynomial fit.
	N_vec : vector_like; Vector which contains the number of fitted data points per source.
	dBIC_Thres : float_like; Delta BIC threshold above which the minimum BIC corresponding to the prefered
	model is significant. If dBIC values of other models are all below this threshold, then the model with
	the least number of parameters is chosen. Typically this value is 2-6. 
	BICcond : Boolean; If False then only the minimum order vector is returned, if True then the BIC_arr is
	also returned.

	# Future versions will incorporate a single source example.
	"""

	# Defining the k*ln(N) array for the BIC:
	k_logN_arr = np.outer(np.log(N_vec),k_vec)

	# Calculating the BIC:
	BIC_arr = Chisqd_arr + k_logN_arr

	# Will index this later to determine the minimum order.
	k_vec_arr = np.ones([len(N_vec),len(k_vec)])*k_vec

	# Finding the minimum BIC value, and the corresponding index vector:
	# Setting all zero values to large numbers for indexing purposes.
	BIC_arr[BIC_arr == 0.0] = 1e+12

	# Vector of the minimum BIC values:
	min_BIC_vec = np.min(BIC_arr, axis=1)
	
	# The delta-BIC array, determined relative to the minimum BIC value:
	dBIC_arr = BIC_arr - min_BIC_vec[:,None]

	# Manipulating the array for indexing purposes:
	dBIC_arr[dBIC_arr >= dBIC_thresh] = dBIC_thresh
	dBIC_arr[np.logical_and(dBIC_arr > 0.0, dBIC_arr < dBIC_thresh)] = 1.0 #This step might not be necessary.

	# Summing the dBIC values:
	dBIC_sum = np.sum(dBIC_arr,axis=1)

	# Where the chosen order is the clearly the best choice:
	min_BIC_True = np.mod(dBIC_sum,dBIC_thresh) == 0.0

	# Where the first model/order is clearly the best choice:
	# We should prefer models with less parameters:
	First_BIC_True = dBIC_arr[:,0] == 0.0

	# The other cases where there is no clear optimal order:
	Other_BIC_True = np.logical_or(First_BIC_True,min_BIC_True) == False

	min_dBIC_ind = np.argmin(dBIC_arr,axis=1)# min dBIC, or where dBIC == 0.

	# In cases where there is no clear significant model, choose the one with less parameters.
	min_dBIC_ind[Other_BIC_True] = min_dBIC_ind[Other_BIC_True] - 1

	# Determining the optimally fit order for each source.
	Order_vec = k_vec[min_dBIC_ind[:,None]] - 1

	if BICcond != False:
		# The case where the user wants to return the BIC array.
		return Order_vec, BIC_arr

	else:
		# The case where the user just wants the order array.
		return Order_vec

if __name__ == "__main__":
	# Plotting stuff:
	#import matplotlib.pyplot as plt
	#import matplotlib.ticker as mtick
	#import matplotlib.gridspec as gridspec
	#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
	#plt.style.use('seaborn-white')
	#plt.rcParams['mathtext.fontset'] = 'stix'
	#plt.rcParams['font.family'] = 'STIXGeneral'
	
	# Defining the parser:
	usage="Usage: %prog [options]\n"
	parser = OptionParser(usage=usage)
	parser.add_option('--freq_cent',dest="freq_cent",default=None,help="Input the normalisation frequency.\n")
	parser.add_option('--catalogue',dest="catalogue",default=None,help="Input list of sources with Zen, Az, RA and Dec.\n")
	parser.add_option('--radec',dest="radec",default=True,help="If False input array only has Zenith and Azimuth columns, else if true the array has both RA, Dec and Zen and Azimuth columns.\n")
	parser.add_option('--suppress_output',dest="Suppress_output",default=False,help="If true then JOOF.py only outputs figures, else JOOF.py prints")

	# Add option for primary beam model inputs.
	(options, args) = parser.parse_args()
	
	# Start of the script:
	start0 = time.time()

	print "\nJOOF.py Version = {0}\n".format(__version__)

	if options.Suppress_output == False:
		print "##############################################################"
		print "# Loading catalogue and initialising parameters:             #"
		print "##############################################################\n"

	# Reading in table of source data/header:
	header = fits.getheader("{0}".format(options.catalogue))
	data = fits.getdata("{0}".format(options.catalogue))

	# Central frequency:
	freq_cent = float(options.freq_cent) # in [MHz]
	
	t = Table(data)#Puts fits data into a readable table format.
	Col_names = t.colnames
	
	if options.Suppress_output == False:
		print "Catalogue = {0}".format(options.catalogue)
		print "Column names = {0}".format(list(Col_names))
	
	# First 4 columns should be Zen, Az, Ra, Dec, the rest of the columns should be
	if options.radec == True:
		chans = Col_names[4:]
	elif options.radec == False:
		chans = Col_names[2:]
	else:
		print "Input error:"
		print "radec should be True or False, default is True."
		print "Catalogue format should be [Zen,Az,RA,DEC,chan_1,...,chan_n]"

	# Converting the channel names to frequency [Hz]:
	Nu_norm = np.array([float(i) for i in chans])*1.28/freq_cent
	
	# Converting table columns into a 2D array:
	Beam_arr = np.array([list(t[chans][j]) for j in range(len(t[chans]))])

	# Specifying the array dimensions:
	k_max = len(chans)-1 # We don't care about the highest order.
	N_sources = len(Beam_arr[:,0]) # Number of sources in the catalogue.
	N_chans = len(Nu_norm) # Number of coarse channels.

	if options.Suppress_output == False:

		print "Number of channels = {0}".format(N_chans)
		print "Number of sources = {0}".format(N_sources)
		print "Max order polynomial to fit = {0}\n".format(k_max-1)
		print "##############################################################"
		print "# Running the Optimal Order Function (OOF):                  #"
		print "##############################################################\n"
		print "Calculating the polynomial coefficients up to order = {0} :\n".format(k_max-1)

	# Determining the polynomial coefficient array:
	#poly_co_arr = poly_coefficients(Nu_norm,np.transpose(Beam_arr),k_max,True)
	#poly_co_arr = poly_coefficients(np.log10(Nu_norm),np.transpose(np.log10(Beam_arr)),k_max)
	poly_co_arr = poly_coefficients(np.log10(Nu_norm),np.transpose(np.log10(Beam_arr)),np.arange(k_max)+1)

	if options.Suppress_output == False:
		print "\nCalculating the polynomial values for N = {0} sources, for N = {1} channels, for each polynomial order.\n".format(N_sources,N_chans)

	# Specifying the polynomial cube:
	poly_ord_arr = np.zeros([k_max,N_sources,N_chans])
	
	# Random sampling for plotting purposes.
	ind = np.random.randint(0,N_sources)
	
	# Calculating the polynomial for source at each order:
	for j in range(k_max):
		if options.Suppress_output == False:
			print "Calculating order {0} polynomial array".format(j)

		poly_ord_arr[j,:,:] = 10**poly_calc(np.log10(Nu_norm),poly_co_arr[j,:,:])
	
	if options.Suppress_output == False:
		print "\nDetermining the optimal order polynomial for each source.\n"
	# Determining the optimal order:
	ord_ar, resid_ar, OOF_ar = OOF(Nu_norm,Beam_arr,poly_ord_arr,figcond=True,\
		Zen=np.array(t['Zen']),theta=np.array(t['theta']),poly_co_list=poly_co_arr)

	end0 = time.time()
	
	if options.Suppress_output == False:
		print "\nTotal runtime = ", np.round(end0-start0,3),"s\n"
		print "##############################################################"
		print "#                            End                             #"
		print "##############################################################"
else:
	# Plotting stuff:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.ticker as mtick
	pass


#def poly_fit_plot(Nu_norm,Beam_vec,poly_co_list,k_vec,y_errors=None,Zen=None,theta=None,radec=False,plot=False):
#	"""
#	This function plots the log-beam response across the bandwidth, and the
#	even or odd (depending on the user input) log-polynomial fits to the 
#	log-beam response. The secondary plot is the squared residuals for each fitted
#	order.
#
#	Args:
#
#	Nu_norm: Vector; Normalised frequecny vector, principally can be any generic x vector.
#	Beam_vec: Vector; Beam value, principally can be any generic y vector.
#	Poly_co_list: Array; Array of polynomial coefficients. Each row is a sucessive order, and each column is 
#	the coefficient for that row in decending order, from the highest to lowerst coefficeint.
#	y_errors: Vector; Vector of associated y-errors.
#	k_vec: Vector; Vector of synonomous to order, whe
#	Zen/Dec: Scalar, default = None; The Azimuth angle. Not None in the non-generic polynomial case. This is
#	the declination when radec = True.	
#	Theta/RA: Scalar, default = None; The Azimuth angle. Not None in the non-generic polynomial case. This is
#	the RA when radec = True.
#	even: Boolean, default = True; Plot only the even polynomials when True,plot odd polynomials when False.
#	radec: Boolean, default = False; Plot local polar coords when False, plot RA/DEC coords when True. 
#	plot: Boolean, default = False; If True, show plot, if False don't show plot.
#	"""
#
#	# Necessary parameters:
#	k_max = len(poly_co_list[:,0])
#
#	# Determine how many parameters there are per entry:
#	poly_k_vec = np.array([len(poly_co_list[i,:][np.abs(poly_co_list[i,:]) > 0.0]) for i in range(len(poly_co_list[:,0]))])
#
#	# Creating an array of indices to plot:
#	order_index_list = np.array([np.arange(len(poly_k_vec))[poly_k_vec == k_vec[i]] for i in range(len(k_vec))]).flatten()
#
#	# Initialising:
#	fig = plt.figure(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
#	frame1 = fig.add_axes((0.1,0.3,0.8,0.6))
#	#frame1.set_xticks([0.95, 1, 1.05])
#	#frame1.tick_params(axis="y", labelsize=15)
#	
#	# Indexing for non-zero beam values:
#	index_list = Beam_vec > 0.0
#
#	# Plotting the initial log-beam response:
#	if np.any(y_errors) != None:
#		# Case where there are input y-errors.
#
#		# Setting errors greater than y to be 0.999*y as to not break the plotting.
#		y_errors_assym = y_errors
#		y_errors_assym[(y_errors-Beam_vec) > 0.0] = 0.999*Beam_vec[(y_errors-Beam_vec) > 0.0]
#
#		plt.errorbar(Nu_norm[index_list],Beam_vec[index_list],yerr=[y_errors_assym[index_list],y_errors[index_list]]\
#			,ecolor='k',fmt='o',markeredgecolor='k',capsize=1.5,elinewidth=0.5,label='chans')
#		# Setting the y-limits
#		plt.ylim([np.min(Beam_vec[index_list])/2.5,2.5*np.max(Beam_vec[index_list])])
#	else:
#		# Case where no errors are input.
#		plt.scatter(Nu_norm[index_list],Beam_vec[index_list],label=r'$\rm{Beam/Chan}$',edgecolors = 'k')
#		# Setting the y-limits
#		plt.ylim([np.min(Beam_vec[index_list])/1.1,1.1*np.max(Beam_vec[index_list])])
#
#	if radec == True and Zen != None and theta != None:
#		# In this case, Zenith and Azimuth are replaced by Dec and Ra:
#		plt.title("RA = {0} deg, DEC = {1} deg".format(np.round(np.degrees(theta),2)\
#			,np.round(np.degrees(Zen),2)),fontsize = 14)
#	elif radec == False and Zen != None and theta != None:
#		# In this case were are dealing with local polar coordinates:
#		plt.title("Azimuth = {0} deg, Zenith = {1} deg".format(np.round(np.degrees(theta),2)\
#			,np.round(np.degrees(Zen),2)),fontsize = 14)
#
#	# Creating an array of smooth points:
#	Nu_norm_smth = np.linspace(np.min(Nu_norm[index_list]),np.max(Nu_norm[index_list]),1000)
#	
#	i = 0
#	for index in order_index_list:
#
#		# Setting the figure:
#		frame1 = fig.add_axes((0.1,0.3,0.8,0.6))
#		plt.loglog()
#		
#		# Fit axis:
#		## Plotting the polynomial fits:
#		plt.plot(Nu_norm_smth,10**poly_calc(np.log10(Nu_norm_smth),poly_co_list[index,:])[0],ls='--',label="order={0}".format(k_vec[i]-1))
#		plt.ylabel(r'$\rm{S_b}$',fontsize = 16)
#		plt.legend()
#		frame1.tick_params(axis="x", labelsize=15)
#		frame1.set_xticks([0.95, 1, 1.05])
#		frame1.tick_params(axis="y", labelsize=15)
#		#frame1.get_xaxis().set_major_formatter(mtick.ScalarFormatter())
#		#frame1.set_xticklabels([])
#	
#		# Residual axis:
#		## Plotting the residual fits:
#		frame2 = fig.add_axes((0.1,0.1,0.8,0.2))
#		plt.loglog()
#		plt.scatter(Nu_norm[index_list],((10**poly_calc(np.log10(Nu_norm)[index_list],poly_co_list[index,:])[0] - Beam_vec[index_list])/Beam_vec[index_list])**2,\
#			edgecolors = 'k',label='order={0}'.format(k_vec[i]-1))
#
#		plt.ylabel(r'$\rm{\chi^2}$',fontsize = 16)
#		plt.xlabel(r'$\rm{\nu/\nu_0}$',fontsize = 16)
#		
#		#frame2.set_xticks([0.95, 1, 1.05])
#		#frame2.tick_params(axis="x", labelsize=15)
#		#frame2.tick_params(axis="y", labelsize=15)
#		
#		#frame2.get_xaxis().set_major_formatter(mtick.ScalarFormatter())
#		#frame2.set_xticklabels([])
#	
#		i += 1 # Counter for indexing the k_vec inside the for loop.
#
#	# The figure is automatically saved:
#	plt.savefig("Poly_fit.png")
#	print "Poly_fit.png saved!"
#
#	if plot == True:
#
#		# Figure is only shown if specified:
#		plt.show()
#	else:
#		plt.close()
#
#	return None
