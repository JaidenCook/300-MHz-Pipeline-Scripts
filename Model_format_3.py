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
matplotlib.use('Agg')
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

# MWA stuff:
from mwa_pb import primary_beam as pb

# Optimal Order Functions:
from JOOF import poly_coefficients, poly_calc, OOF

def mwa_alt_az_za(obsid, ra=None, dec=None, degrees=False):
	"""
	Calculate the altitude, azumith and zenith for an obsid
	
	Args:
	obsid : The MWA observation id (GPS time)
	ra : The right acension in HH:MM:SS
	dec : The declintation in HH:MM:SS
	degrees: If true the ra and dec is given in degrees (Default:False)
	"""
	from astropy.time import Time
	from astropy.coordinates import SkyCoord, AltAz, EarthLocation
	from astropy import units as u
	
	obstime = Time(float(obsid),format='gps')
	
	if ra is None or dec is None:
	#if no ra and dec given use obsid ra and dec
		ra, dec = get_common_obs_metadata(obsid)[1:3]
	
	if degrees:
		sky_posn = SkyCoord(ra, dec, unit=(u.deg,u.deg))
	else:
		sky_posn = SkyCoord(ra, dec, unit=(u.hourangle,u.deg))
	earth_location = EarthLocation.of_site('Murchison Widefield Array')
	#earth_location = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
	altaz = sky_posn.transform_to(AltAz(obstime=obstime, location=earth_location))
	Alt = altaz.alt.deg
	Az = altaz.az.deg
	Za = 90. - Alt
	return Alt, Az, Za


if __name__ == "__main__":
	# Timing the script duration.
	start0 = time.time()
	
	# Parser options:
	from optparse import OptionParser
	
	usage="Usage: %prog [options]\n"
	parser = OptionParser(usage=usage)
	parser.add_option('--obsid',dest="obsid",default=None,help="Input OBSID")
	parser.add_option('--freq',dest="freq",default=None,help="Input frequency")
	parser.add_option('--delays',dest="delays",default=None,help="Input delays")
	parser.add_option('--catalogue',dest="catalogue",default=None,help="Input GLEAM catalogue to use")
	parser.add_option('--ra',dest="ra_pnt",default=None,help="RA of pointing centre (deg)")
	parser.add_option('--dec',dest="dec_pnt",default=None,help="Dec of pointing centre (deg)")
	parser.add_option('--threshold',dest="threshold",default=0.0127,help="Input threshold sigma scalar")
	
	# Add option for primary beam model inputs -- make default a typical GLEAMY2 repository location
	(options, args) = parser.parse_args()
	
	# Loading in parameters from the metafits file:
	metafits = "{0}.metafits".format(int(options.obsid))
	hdus = fits.open(metafits)
	meta_dat = hdus[0].header
	
	# Reading inputs from the command line.
	Tot_Sky_Mod = options.catalogue
	
	# Reading in the delay string, converting to a list where each element is a float.
	delays = options.delays
	delays = delays[1:len(delays)-1].split(",")
	delays_flt = [float(i) for i in delays]
	print "Delays = {0}".format(delays_flt)
	
	# Initialising other parameters.
	obsid = float(options.obsid)
	cent_freq = (float(options.freq))*1e+6 #Central frequency
	ra_pnt = float(options.ra_pnt)
	dec_pnt = float(options.dec_pnt)
	threshold = float(options.threshold)
	
	if not os.path.exists(Tot_Sky_Mod):
	   print "Can't find the total sky model catalogue."
	   sys.exit(1)
	
	freq_str = "%03.0f" % (cent_freq/1e6)
	print "Central frequency = {0} MHz".format(freq_str)
	
	# Loading in the channels for a given observation.
	chans = meta_dat['CHANNELS']
	print "Channels = {0}\n".format(chans)
	chans = chans[:len(chans)].split(",")
	Nu = (np.array([float(i) for i in chans])*1.28)*1e+6 #Converting to hz
		
	# Add columns with RA and Dec in sexagesimal format
	os.system('stilts tpipe in='+Tot_Sky_Mod+' cmd=\'addcol ra_str "degreesToHms(RA,2)"\' cmd=\'addcol dec_str "degreesToDms(DEC,2)"\' out=tmp.fits')
	
	# Read in GLEAM IDR4
	temp_table = fits.open('tmp.fits')[1].data
	
	RA_samp = temp_table['RA']
	DEC_samp = temp_table['DEC']
	
	#"""
	################################################################################
	# Thresholding all sources above the horizon, converting RA/DEC to Az/Alt
	################################################################################
	#"""
	
	# Determining the Altitude, Azimuth and Zenith for each source in the total catalogue.
	Alt0, Az0, Zen0 =  mwa_alt_az_za(obsid, RA_samp, DEC_samp, True)
	
	Alt0_samp = Alt0[Alt0 > 0.0]
	Az0_samp = np.radians(Az0[Alt0 > 0.0])
	Zen0_samp = (pi/2) - np.radians(Alt0_samp)
	
	# Trying to solve an indexing issue with mwa_pb.
	Az0_samp = [Az0_samp]
	Zen0_samp = [Zen0_samp]
	
	print "Number of sources above the horizon:",len(Alt0_samp)
	
	# Creating new sampled table.
	Sky_Mod = temp_table[Alt0 > 0.0]
	
	#"""
	################################################################################
	# Thresholding the apparent brightest 1500 sources in the OBSID
	################################################################################
	#"""
	
	# Determining the beam power at the central frequency for each source above the horizon.
	beam_power = pb.MWA_Tile_full_EE(Zen0_samp,Az0_samp,cent_freq,delays_flt)
	
	# Splitting the xx and yy cross correlation power components
	beam_power_XX = np.array(beam_power[0])
	beam_power_YY = np.array(beam_power[1])
	
	# Determining the average power or stokes I power.
	beamvalue = (beam_power_XX + beam_power_YY)/2
	
	S300_fit = Sky_Mod.field("Fint300")
	
	# Determining the apparent flux density for each source.
	S_centralfreq = S300_fit
	S_centralfreq_uncorrected = S_centralfreq*beamvalue[0,:]
	
	# Thresh_indices = S_centralfreq_uncorrected >= 5*threshold
	# This will make calibration go faster.
	Thresh_indices = np.argsort(S_centralfreq_uncorrected)[len(S_centralfreq_uncorrected)-1500:]
	
	# Subsetting by applying flux threshold cut:
	S_centralfreq = S_centralfreq[Thresh_indices]
	beamvalue = beamvalue[0,:][Thresh_indices]
	S300_fit = S300_fit[Thresh_indices]
	
	print "Number of sources after thresholding = # ", len(S_centralfreq_uncorrected)
	
	# Loading in the other vectors:
	name = Sky_Mod.field("Name")[Thresh_indices]
	ra = Sky_Mod.field("RA")[Thresh_indices]
	dec = Sky_Mod.field("DEC")[Thresh_indices]
	ra_str = Sky_Mod.field("ra_str")[Thresh_indices]
	dec_str = Sky_Mod.field("dec_str")[Thresh_indices]
	int_flux_wide = Sky_Mod.field("Fint300")[Thresh_indices]
	a_wide = Sky_Mod.field("Major")[Thresh_indices]
	b_wide = Sky_Mod.field("Minor")[Thresh_indices]
	pa_wide = Sky_Mod.field("PA")[Thresh_indices]
	peak_flux_wide = Sky_Mod.field("Fint300")[Thresh_indices]
	flags = Sky_Mod.field("flag")[Thresh_indices]
	coeff_arr = np.array(Sky_Mod.field("coefficients"))[Thresh_indices,:]
	os.remove('tmp.fits')
		
	#"""
	################################################################################
	# Fitting the beam curvature
	################################################################################
	#"""
	
	# Initialising beam_cube:
	flat_beam_cube = np.empty([len(Nu),len(name)])
	
	# Retrieving the Az and Zen angles for the top 1500 sources.
	Az0_samp = [np.array(Az0_samp)[0,:][Thresh_indices]]
	Zen0_samp = [np.array(Zen0_samp)[0,:][Thresh_indices]]
	
	# Calculating the beam for each source, and each coarse channel.
	for i in range(len(Nu)):
		print "Generating channel {0} beam".format(chans[i])
		temp_beam_power = pb.MWA_Tile_full_EE(Zen0_samp,Az0_samp,Nu[i],delays_flt)
		flat_beam_cube[i,:] = ((np.array(temp_beam_power[0]) + np.array(temp_beam_power[1]))/2.0).flatten()
	
	#"""
	################################################################################
	# Determining the Polynomial Coefficeints:
	################################################################################
	
	print "Determining the beam spectral index and curvature"
	
	# Important parameters:
	k_max = 12
	N_sources = len(name)
	N_chans = len(Nu)
	
	# Determining the polynomial coefficient array.
	poly_co_arr = poly_coefficients(np.log10(Nu/cent_freq),np.log10(flat_beam_cube),np.arange(k_max)+1)
	
	# Specifying the polynomial cube:
	poly_ord_arr = np.zeros([k_max,N_sources,N_chans])
	
	# Generating random integer for plotting purposes:
	ind = np.random.randint(0,N_sources)
	
	# Determining the polynomial values for each order, for each coarse channel for each source.
	for j in range(k_max):
		print "Calculating order {0} polynomial array".format(j)
		poly_ord_arr[j,:,:] = 10**poly_calc(np.log10(Nu/cent_freq),poly_co_arr[j,:,:])
	
	# Determining the optimal order:
	ord_ar, resid_ar, OOF_ar = OOF(Nu/cent_freq,np.transpose(flat_beam_cube),poly_ord_arr,\
		figcond=True,Zen=Zen0_samp[0],theta=Az0_samp[0],poly_co_list=poly_co_arr)
	
	# Getting the coefficients corresponding to the optimal order:
	temp_poly_co = poly_co_arr[ord_ar,np.arange(N_sources),:]
	
	# At the central frequency all the log polynomial terms are 0.
	beamvalue_approx = 10**temp_poly_co[:,k_max - 1]
	
	# Apparent flux:
	S_centralfreq_uncorrected = S_centralfreq*beamvalue_approx
	
	# Specifying source polynomial coefficients, without the log-central fllux density:
	source_poly_co_ar = np.zeros([N_sources,k_max-1])
	
	# Apparent polynomial coefficients:
	app_poly_co_ar = np.zeros(np.shape(temp_poly_co))# Initialising.
	
	app_poly_co_ar[:,len(coeff_arr[0,:])-k_max:] = np.around((temp_poly_co[:,len(coeff_arr[0,:])-k_max:] + coeff_arr[:,::-1]),decimals=2)
	
	# Converting columns in the coefficient array into tuples, and concatenating them into a list.
	app_poly_co_tup_list = [tuple(app_poly_co_ar[l,:]) for l in range(N_sources)]
	
	################################################################################
	# Creating output VO table
	################################################################################
	#"""
	
	print "Writing V0 table"
	
	filename_body='model'
	
	newvot = Table( [name, ra_str, dec_str, ra, dec, coeff_arr[:,1], S_centralfreq, beamvalue, beamvalue_approx, S_centralfreq_uncorrected,\
		app_poly_co_tup_list,peak_flux_wide, int_flux_wide, a_wide, b_wide, pa_wide],\
	    names=('Name','ra_str','dec_str','ra','dec', 'alpha','S_centralfreq',
	        'beamvalue','beamvalue_approx','S_centralfreq_uncorrected','apparent_poly_coeff','peak_flux_wide','int_flux_wide',\
	        'a_wide', 'b_wide', 'pa_wide'), meta={'name': 'first table'} )
	
	writetoVO(newvot, filename_body+"_morecolumns_temp.vot")
	
	#"""
	print 'Completed'
	end0 = time.time()
	
	print "Total runtime = ", np.round(end0 - start0,3),"s\n"
