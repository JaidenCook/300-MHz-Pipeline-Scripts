import os
import sys
import time
from datetime import datetime
import glob
import shutil
import re
from math import pi

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

import numpy as np

from scipy.optimize import leastsq
from scipy.optimize import curve_fit

# Astropy stuff
from astropy import wcs
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import Table,Column
from astropy.io.votable import writeto as writetoVO

from mwa_pb import primary_beam as pb

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

def fit(nu,Stemp,init,mod,Serrtemp=None):
	"""
	Defining a fitting function for parallel processing.
	"""
	#log_Serrtemp = 1

	if Serrtemp is None:
		log_Serrtemp = 1
	else:
		log_Serrtemp = np.log(10)*Serrtemp/Stemp
		
	c1_err = []
	out = leastsq(mod,init,args=(nu,np.log10(Stemp),log_Serrtemp),full_output=1, epsfcn=0.0001)
	c1 = out[0]
	cov = out[1]

	for l in range(len(c1)):
		c1_err.append(np.sqrt(cov[l,l]))
	
	return c1

################################################################################
#Setting up the parser
################################################################################
from optparse import OptionParser

usage="Usage: %prog [options]\n"
parser = OptionParser(usage=usage)
parser.add_option('--obsid',dest="obsid",default=None,help="Input OBSID")
parser.add_option('--size',dest="size",default=None,help="Input image dimensions")
parser.add_option('--reztable',dest="reztable",default=True,help="True/False, output high residual fringe table")

(options, args) = parser.parse_args()
################################################################################
#Initialsing
################################################################################
#
#Initialising output file.
newfile = open("{0}_beam_residual_o.txt".format(float(options.obsid)),"w")

#Lambda functions used for fitting the beam:
log_Quart = lambda a,nu: a[0] + a[1]*(np.log10(nu/300e+6)) + a[2]*(np.log10(nu/300e+6))**2 + a[3]*(np.log10(nu/300e+6))**3 + a[4]*(np.log10(nu/300e+6))**4
log_Cub = lambda a,nu: a[0] + a[1]*(np.log10(nu/300e+6)) + a[2]*(np.log10(nu/300e+6))**2 + a[3]*(np.log10(nu/300e+6))**3
log_Quad = lambda a,nu: a[0] + a[1]*(np.log10(nu/300e+6)) + a[2]*(np.log10(nu/300e+6))**2
log_Lin = lambda a,nu: a[0] + a[1]*(np.log10(nu/300e+6))
Sqd_Quart_Residuals = lambda a,nu,Y,Yerr: ((Y-log_Quart(a,nu))**2)/Yerr**2
Sqd_Cub_Residuals = lambda a,nu,Y,Yerr: ((Y-log_Cub(a,nu))**2)/Yerr**2
Sqd_Quad_Residuals = lambda a,nu,Y,Yerr: ((Y-log_Quad(a,nu))**2)/Yerr**2
Sqd_Lin_Residuals = lambda a,nu,Y,Yerr: ((Y-log_Lin(a,nu))**2)/Yerr**2

# Loading in parameters from the metafits file:
# These will become parsed parameters.
metafits = "{0}.metafits".format(int(options.obsid))
hdus = fits.open(metafits)
meta_dat = hdus[0].header
img_dim = int(options.size)

print "OBSID = {0}".format(int(options.obsid))
print "Image dimensions = {0}x{0}".format(img_dim)
newfile.write("OBSID = {0}\n".format(int(options.obsid)))
newfile.write("Image dimensions = {0}x{0}\n".format(img_dim))

# Loading in the channels for a given observation.
chans = meta_dat['CHANNELS']
print "Channels = {0}\n".format(chans)
newfile.write("Channels = {0}\n".format(chans))
chans = chans[:len(chans)].split(",")
Nu = (np.array([float(i) for i in chans])*1.28)*1e+6 #Converting to hz

# Loading in the delays for that given observation pointing.
delays = meta_dat['DELAYS']
print "Delays = {0}".format(delays)
newfile.write("Delays = {0}\n".format(delays))
delays = delays[:len(delays)].split(",")
delays_flt = np.array([float(j) for j in delays])

# Loading in the central frequency of the observation.
cent_freq = meta_dat['FREQCENT']*1e+6
print "Central frequency = {0} Hz".format(cent_freq)
newfile.write("Central frequency = {0} Hz\n".format(cent_freq))

# Initialising beam_cube:
flat_beam_cube = np.empty([len(Nu),img_dim*img_dim])

# Initialising angular coordinates.
Alt0_samp = np.linspace(0,pi/2,img_dim)
Az0_samp = np.linspace(0,2*pi,img_dim)
Zen0_samp = (pi/2) - Alt0_samp

# Creating a grid of values:
Zen,th = np.meshgrid(Zen0_samp, Az0_samp)

for i in range(len(Nu)):
	print "Generating channel {0} beam".format(chans[i])
	temp_beam_power = pb.MWA_Tile_full_EE(Zen,th,Nu[i],delays_flt)
	flat_beam_cube[i,:] = ((np.array(temp_beam_power[0]) + np.array(temp_beam_power[1]))/2.0).flatten()

################################################################################
#Fitting log-polynomial models
################################################################################
#Quartic model fit.
print "Fitting quartic model"
start = time.time()
results_quart = np.array(Parallel(n_jobs=-1)(delayed(fit)(Nu,flat_beam_cube[:,h],[1,1,1,1,1],Sqd_Quart_Residuals) for h in tqdm(range(img_dim*img_dim))))
end = time.time()

print "Quartic Fitting Duration = ", np.round(end - start,3),"s\n"
newfile.write("Quartic fit duration = {0} s\n".format(np.round(end - start,3)))
print "Fitting cubic model"

#Cub model fit.
start = time.time()
results_cub = np.array(Parallel(n_jobs=-1)(delayed(fit)(Nu,flat_beam_cube[:,i],[1,1,1,1],Sqd_Cub_Residuals) for i in tqdm(range(img_dim*img_dim))))
end = time.time()

print "Cubic Fitting Duration = ", np.round(end - start,3),"s\n"
newfile.write("Cubic fit duration = {0} s\n".format(np.round(end - start,3)))
print "Fitting quadratic model"

#Quadratic model fit.
start = time.time()
results_quad = np.array(Parallel(n_jobs=-1)(delayed(fit)(Nu,flat_beam_cube[:,j],[1,1,1],Sqd_Quad_Residuals) for j in tqdm(range(img_dim*img_dim))))
end = time.time()

print "Quadratic Fitting Duration = ", np.round(end - start,3),"s\n"
newfile.write("Quadratic fit duration = {0} s\n".format(np.round(end - start,3)))
print "Fitting linear model"

#Linear model parallel fit, used for comparison.
start = time.time()
results_lin = np.array(Parallel(n_jobs=-1)(delayed(fit)(Nu,flat_beam_cube[:,k],[1,1],Sqd_Lin_Residuals) for k in tqdm(range(img_dim*img_dim))))
end = time.time()

print "Linear Fitting Duration = ", np.round(end - start,3), "s"
newfile.write("Linear fit duration = {0} s\n".format(np.round(end - start,3)))

# Setting the coefficient maps:
Slog_lin = results_lin[:,0].reshape(img_dim,img_dim)
alpha_lin = results_lin[:,1].reshape(img_dim,img_dim)

Slog_quad = results_quad[:,0].reshape(img_dim,img_dim)
alpha_quad = results_quad[:,1].reshape(img_dim,img_dim)
q_quad = results_quad[:,2].reshape(img_dim,img_dim)

Slog_cub = results_cub[:,0].reshape(img_dim,img_dim)
alpha_cub = results_cub[:,1].reshape(img_dim,img_dim)
q_cub = results_cub[:,2].reshape(img_dim,img_dim)
r_cub = results_cub[:,3].reshape(img_dim,img_dim)

Slog_quart = results_quart[:,0].reshape(img_dim,img_dim)
alpha_quart = results_quart[:,1].reshape(img_dim,img_dim)
q_quart = results_quart[:,2].reshape(img_dim,img_dim)
r_quart = results_quart[:,3].reshape(img_dim,img_dim)
s_quart = results_quart[:,4].reshape(img_dim,img_dim)

# Creating beam models from the determined beam coefficients:
S300quart_mod = 10**log_Lin([Slog_quart,alpha_quart,q_quart,r_quart,s_quart],cent_freq) 
S300cub_mod = 10**log_Lin([Slog_cub,alpha_cub,q_cub,r_cub],cent_freq)
S300quad_mod = 10**log_Quad([Slog_quad,alpha_quad,q_quad],cent_freq)
S300lin_mod = 10**log_Lin([Slog_lin,alpha_lin],cent_freq)

# Creating list of coefficient maps:
Coeff_list = [Slog_lin,alpha_lin,Slog_quad,alpha_quad,q_quad,Slog_cub,alpha_cub,q_cub,\
r_cub,Slog_quart,alpha_quart,q_quart,r_quart,s_quart]

# Creating 300MHz beam model:
temp_beam_power = pb.MWA_Tile_full_EE(Zen,th,cent_freq,delays_flt)
S300_beam_mod = ((np.array(temp_beam_power[0]) + np.array(temp_beam_power[1]))/2.0)

# Calculating the residual maps:
S300_quart_res =  (S300_beam_mod - S300quart_mod)/S300_beam_mod
S300_cub_res =  (S300_beam_mod - S300cub_mod)/S300_beam_mod
S300_quad_res =  (S300_beam_mod - S300quad_mod)/S300_beam_mod
S300_lin_res = (S300_beam_mod - S300lin_mod)/S300_beam_mod


################################################################################
#Creating Fringe table:
################################################################################
if options.reztable  == 'True':
	Fringe_ind = np.array(np.abs(S300_quart_res) > 0.6).flatten()
	
	S_test = S300_quart_res[np.abs(S300_quart_res) > 0.6]
	print "Number of high residual fringe sources = ", len(S_test)
	
	Fringe_list = [Zen.flatten()[Fringe_ind],th.flatten()[Fringe_ind]]
	for i in range(len(chans)):
		Fringe_list.append(flat_beam_cube[i,Fringe_ind])
	
	for j in range(len(Coeff_list)):
		Fringe_list.append(Coeff_list[j].flatten()[Fringe_ind])
	
	Col_names = ['Zen','theta']
	Coeff_names = ['S1','a1','S2','a2','q2','S3','a3','q3','r3','S4','a4','q4','r4','u4']
	
	# Creating a string of column names:
	Col_names = Col_names + chans + Coeff_names
	
	print "Creating output high residual fringe table!"
	Output_Table = Table(Fringe_list,names=Col_names,meta={'name':'first table'})
	Output_Table.write("{0}_fringe_table.fits".format(int(options.obsid)),"w")
else:
	print "Output table is False"

################################################################################
#Plotting residual maps
################################################################################
print "Plotting residual maps"

r = np.sin(Zen)
fig = plt.figure(figsize = (12.5,9.5), dpi=90)

ax1 = fig.add_subplot(221,projection="polar")
pcm1 = ax1.pcolormesh(th, r, S300_lin_res, cmap='jet')
ax1.title.set_text("Linear Residuals")
ax1.title.set_position([0.5, 1.05])
ax1.set_yticklabels([])
ax1.set_theta_zero_location('N')
cbaxes1 = fig.add_axes([0.465, 0.55, 0.01, 0.38]) 
fig.colorbar(pcm1, ax=ax1, cax=cbaxes1)


ax2 = fig.add_subplot(222,projection="polar")
pcm2 = ax2.pcolormesh(th, r, S300_quad_res, cmap='jet')
ax2.title.set_text("Quadratic Residuals")
ax2.title.set_position([0.5, 1.05])
ax2.set_yticklabels([])
ax2.set_theta_zero_location('N')
cbaxes2 = fig.add_axes([0.935, 0.55, 0.01, 0.38]) 
fig.colorbar(pcm2, ax=ax2, cax=cbaxes2)

ax3 = fig.add_subplot(223,projection="polar")
pcm3 = ax3.pcolormesh(th, r, S300_cub_res, cmap='jet')
ax3.title.set_text("Cubic Residuals")
ax3.title.set_position([0.5, 1.06])
ax3.set_yticklabels([])
ax3.set_theta_zero_location('N')
cbaxes3 = fig.add_axes([0.465, 0.05, 0.01, 0.38]) 
fig.colorbar(pcm3, ax=ax3, cax=cbaxes3)

ax4 = fig.add_subplot(224,projection="polar")
pcm4 = ax4.pcolormesh(th, r, S300_quart_res, cmap='jet')
ax4.title.set_text("Quartic Residuals")
ax4.title.set_position([0.5, 1.06])
ax4.set_yticklabels([])
ax4.set_theta_zero_location('N')
cbaxes4 = fig.add_axes([0.935, 0.05, 0.01, 0.38]) 
fig.colorbar(pcm4, ax=ax4, cax=cbaxes4)

plt.tight_layout()
plt.savefig("Residuals_map.png")
#plt.show()

print "Mean quartic residuals = ",np.round(np.mean(np.abs(S300_quart_res)),4)
print "Mean cubic residuals = ",np.round(np.mean(np.abs(S300_cub_res)),4)
print "Mean quadratic residuals = ",np.round(np.mean(np.abs(S300_quad_res)),4)
print "Mean linear residuals = ",np.round(np.mean(np.abs(S300_lin_res)),4)

print "Std quartic residuals = ",np.round(np.std(np.abs(S300_quart_res)),4)
print "Std cubic residuals = ",np.round(np.std(np.abs(S300_cub_res)),4)
print "Std quadratic residuals = ",np.round(np.std(np.abs(S300_quad_res)),4)
print "Std linear residuals = ",np.round(np.std(np.abs(S300_lin_res)),4)


#Writing important output to file:
newfile.write("Mean quartic residuals = {0}\n".format(np.round(np.mean(np.abs(S300_quart_res)),4)))
newfile.write("Mean cubic residuals = {0}\n".format(np.round(np.mean(np.abs(S300_cub_res)),4)))
newfile.write("Mean quadratic residuals = {0}\n".format(np.round(np.mean(np.abs(S300_quad_res)),4)))
newfile.write("Mean linear residuals = {0}\n".format(np.round(np.mean(np.abs(S300_lin_res)),4)))

newfile.write("Std quartic residuals = {0}\n".format(np.round(np.std(np.abs(S300_quart_res)),4)))
newfile.write("Std cubic residuals = {0}\n".format(np.round(np.std(np.abs(S300_cub_res)),4)))
newfile.write("Std quadratic residuals = {0}\n".format(np.round(np.std(np.abs(S300_quad_res)),4)))
newfile.write("Std linear residuals = {0}\n".format(np.round(np.std(np.abs(S300_lin_res)),4)))

newfile.close()
