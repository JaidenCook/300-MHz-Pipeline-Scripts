#! /usr/bin/env python
#!/pawsey/cle52up04/apps/gcc/4.3.4/python/2.6.9/bin/python

# Take a catalogue (e.g. PUMA crossmatched with GLEAM), v2 - that is in .fits format
# v2 - Reading in more columns so that output can be run through vo2model.py
# Find the sources within a radius of RA and Dec of an image
# Attenuate the sources by the primary beam model
# v1 - Write the sources that have been selected (i.e. those within radius and above the threshold) to a VO table
# [Put ratios down on a grid based on the image WCS - not actually doing this]
version='v6'

# Email sarah.white@icrar.org if stuck!

import os
import sys
import glob
import shutil
import re
import time
import numpy as np
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

from astropy import wcs
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import Table,Column
from astropy.io.votable import writeto as writetoVO

import matplotlib as mpl
mpl.use('Agg') # So does not use display
import matplotlib.pylab as plt

### Beam_Norm() was written by Jaiden Cook
def Beam_Norm(Beam):
    """
    This function renormalises the tile beam model for an arbitrary subband
    from zenith to the maximum sensitivty.
    """
    #Currently set to not normalise the beam.
    FITSData = Beam[0].data
    #tmp.info()
    #Cuts the data from the FITS file into a 2D array.
    #Map = FITSData[0,0,:,:]
    
    #Defining a mask array where the Nan values are true and the non nan values are false.
    #Nan_arr = np.isnan(Map)
    #Need to convert nan values to 0 in order to find the max pixel value.
    #Map = np.nan_to_num(Map)
    #Scale down the beam relative to the maximum.
    #Map = Map/np.max(Map)
    #print "Maximum beam value = ", np.max(Map)
    #Use Nan_arr as a mask to replace the values with numpy nan values.
    #np.place(Map,Nan_arr,np.nan)
    #Rewriting the scaled data to the fits data, then write a new fits file.
    #FITSData[0,0,:,:] = Map
    
    return FITSData

### Fiterr() was written by Jaiden Cook
def Fiterr(fitfunc,init,X,Y,Yerr):
	"""
	This function is designed ot re-run fitting processes for a number of iterations until a convergent answer is achieved. This is not necessarily necessary for all fits, but has been found to give more accurate parameter values when using previous fit params as guesses.
	"""
	c1_err = []
	for j in range(10):
		if j == 0:
			out = leastsq(fitfunc,init,args=(X,Y,Yerr),full_output=1, epsfcn=0.0001)
			c1 = out[0]
		else:
			out = leastsq(fitfunc,init,args=(X,Y,Yerr),full_output=1, epsfcn=0.0001)
			c1 = out[0]
			#cov = out[1]

	#for l in range(len(c1)):
		#c1_err.append(np.sqrt(cov[l,l]))
	
	return c1#, c1_err

### gcd() was written by Paul Hancock, modified by Jaiden Cook
# The following functions are explained at http://www.movable-type.co.uk/scripts/latlong.html
# JC: I vectorised the function to eliminate unnecessary for loops.
def gcd(ra1, dec1, ra2, dec2):
    """
    Great circle distance as calculated by the haversine formula
    ra/dec in degrees
    returns:
    sep in degrees
    """
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(np.radians(dlat)/2)**2
    a += np.cos(np.radians(dec1))*np.cos(np.radians(dec2))*\
np.sin(np.radians(dlon)/2)**2

    a_sqrt = np.sqrt(a)
    a_sqrt[a_sqrt > 1.0] = 1.0 
    sep = np.degrees(2*np.arcsin(a_sqrt))
    return sep

start0 = time.time()

#Lambda functions used for fitting the beam at 200/300MHz:
#
# This will have to be generalised by changing nu= nu/nu_0.
#log_Quad = lambda a,nu: a[0] + a[1]*(np.log10(nu/200e+6)) + a[2]*(np.log10(nu/200e+6))**2
# x = nu/nu_0, generalised for GLEAM bands.
log_Quad = lambda a,x: a[0] + a[1]*(np.log10(x)) + a[2]*(np.log10(x))**2
#Sqd_Quad_Residuals = lambda a,nu,Y,Yerr: ((Y-log_Quad(a,nu))**2)/Yerr**2
Sqd_Quad_Residuals = lambda a,x,Y,Yerr: ((Y-log_Quad(a,x))**2)/Yerr**2

from optparse import OptionParser

usage="Usage: %prog [options]\n"
parser = OptionParser(usage=usage)
parser.add_option('--pbeam',dest="pbeam",default=None,help="Input Stokes I primary beam for 30.72 MHz band")
parser.add_option('--catalogue',dest="catalogue",default=None,help="Input GLEAM catalogue to use")
parser.add_option('--ra',dest="ra_pnt",default=None,help="RA of pointing centre (deg)")
parser.add_option('--dec',dest="dec_pnt",default=None,help="Dec of pointing centre (deg)")
parser.add_option('--threshold',dest="threshold",default=None,help="Input threshold sigma scalar")
# Add option for primary beam model inputs -- make default a typical GLEAMY2 repository location
(options, args) = parser.parse_args()

# This section is required for fitting a quadratic to the beam.
obsid = options.pbeam[0:10]

# The subband beams should be created using build_model.sh.
Beams = ["{0}-0000-psf_beamI.fits".format(obsid),"{0}-0001-psf_beamI.fits".format(obsid),"{0}-0002-psf_beamI.fits".format(obsid),"{0}-0003-psf_beamI.fits".format(obsid)]

# Double checking the Obsid variable is what we exepect.
print "Obsid: ",obsid

# Loading in the central frequencies for each subband.
Nu = np.zeros(len(Beams))
Beam_Cube = []
for i in range(len(Beams)):
	Nu[i]=fits.getheader(Beams[i])['CRVAL3']
	Beam = fits.open(Beams[i])
	Beam_Cube.append(Beam_Norm(Beam)[0,0,:,:])

# This will be used to fit the log-quadratic coefficients.
Beam_Cube = np.array(Beam_Cube)

# Reading inputs from the command line.
beamimage=options.pbeam
GLEAM=options.catalogue
ra_pnt=float(options.ra_pnt)
dec_pnt=float(options.dec_pnt)
threshold=float(options.threshold)
if not os.path.exists(beamimage):
   print "Can't find the specified primary-beam model for 30.72 MHz band."
   sys.exit(1)
if not os.path.exists(GLEAM):
   print "Can't find the GLEAM catalogue."
   sys.exit(1)

# You will have defined beamimage earlier
print "beamimage = ",beamimage
hdulist = fits.open(beamimage)
# Here the beam is normalised from zenith to the maximum sensitivity of the beam.
hdulist[0].data = Beam_Norm(hdulist)

# Get WCS for beam image
w = wcs.WCS(hdulist[0].header,naxis=2) 

# Get image dimensions
xdim = fits.getheader(beamimage)['NAXIS1']
ydim = fits.getheader(beamimage)['NAXIS2']
print 'xdim =',xdim
print 'ydim =',ydim

# Get cell size, in deg (assume square)
cell = fits.getheader(beamimage)['CDELT2']
print 'Cell size (deg) =',cell

# Set radius to half the image size
radius=cell*ydim/2.0 #Fixing to width of 300 MHz primary beam.
#radius=13.6 # added extra 2 degrees to account for off centre pb normalisation.
#radius = 110
print 'Radius (deg) =',radius

# Get central frequency of the various sub-bands
freq = fits.getheader(beamimage)['CRVAL3']   # Units of Hz
print "Frequency test print = ",freq

freq_str = "%03.0f" % (freq/1e6)

# Select sub-band flux which is closest in frequency to the central frequency of the image
freq_list=np.array([76.0,84.0,92.0,99.0,107.0,115.0,122.0,130.0,143.0,151.0,\
158.0,166.0,174.0,181.0,189.0,197.0,204.0,212.0,220.0,227.0,288.0,296.0,300.0,\
304.0,312.0])

freq_list2=['076','084','092','099','107','115','122','130','143','151','158',\
'166','174','181','189','197','204','212','220','227','288','296','300','304',\
'312']

subband_freq=freq_list[np.argmin(abs(freq_list-freq/1e6))]
subband_freq_str=freq_list2[np.argmin(abs(freq_list-freq/1e6))]
subband_flux='Fint'+str(subband_freq_str)
print 'subband_freq:',subband_freq
print 'Sub-band flux which is closest in frequency to central frequency of image:',subband_flux

# Add columns with RA and Dec in sexagesimal format
os.system('stilts tpipe in='+GLEAM+' cmd=\'addcol ra_str "degreesToHms(RAJ2000,2)"\' cmd=\'addcol dec_str "degreesToDms(DEJ2000,2)"\' out=tmp.fits')

# Read in GLEAM IDR4
GLEAMtable = fits.open('tmp.fits')[1].data

name=GLEAMtable.field("Name")
ra=GLEAMtable.field("RAJ2000")
dec=GLEAMtable.field("DEJ2000")
ra_str=GLEAMtable.field("ra_str")
dec_str=GLEAMtable.field("dec_str")
S200_fit=GLEAMtable.field("S_200")
alpha=GLEAMtable.field("alpha")
q_curve=GLEAMtable.field("beta")
#int_flux_subband=GLEAMtable.field(subband_flux)
int_flux_wide=GLEAMtable.field("S_200")
a_wide=GLEAMtable.field("a")
b_wide=GLEAMtable.field("b")
pa_wide=GLEAMtable.field("pa")
peak_flux_wide=GLEAMtable.field("S_200")
os.remove('tmp.fits')

start = time.time()
ang_sep = gcd(ra_pnt,dec_pnt,ra,dec)
end = time.time()

indices=np.where(ang_sep < radius)

################################################################################
#Radial Filtering
################################################################################
name_withinradius=name[indices]
ra_withinradius=ra[indices]
dec_withinradius=dec[indices]
ra_str_withinradius=ra_str[indices]
dec_str_withinradius=dec_str[indices]
S200_fit_withinradius=S200_fit[indices]
alpha_withinradius=alpha[indices]
q_curve_withinradius=q_curve[indices]
#int_flux_subband_withinradius=int_flux_subband[indices]
int_flux_wide_withinradius=int_flux_wide[indices]
a_wide_withinradius=a_wide[indices]
b_wide_withinradius=b_wide[indices]
pa_wide_withinradius=pa_wide[indices]
peak_flux_wide_withinradius=peak_flux_wide[indices]

print "Number of sources after radial filtering = ",len(name_withinradius)

# Next step is to calculate the S_central frequency:
#
# Determining the log-quadratic coefficients for each source:
C1 = np.log10(S200_fit_withinradius)
C2 = alpha_withinradius
C3 = q_curve_withinradius

# This process has been vectorised for efficiency:
S_centralfreq_withinradius = 10**(C1 + C2*(np.log10(freq/200e+6)) + C3*(np.log10(freq/200e+6))**2)

alpha_withinradius = C2 + 2*C3*np.log10(freq/200e+6)

# Now that you have a list of RAs and Decs, find the pixel co-ordinates
imagedata=hdulist[0].data

# This position *should* be at the centre of the image
ra_test=ra_pnt
dec_test=dec_pnt
px, py = w.wcs_world2pix(ra_test, dec_test,1) # The WCS details are already 'contained' within the 'w'. Just need to specify the RA and Dec (in degrees I assume)

print 'Pixel corresponding to the reference position:'
print 'px = ',px,'  py = ',py  # Note that these are decimals, so would need to be rounded before being used as pixel indices

# Check that I get the RA/Dec back if I go the other way round
ra_out, dec_out = w.wcs_pix2world(px, py, 1)

print 'The following should match ra_pnt and dec_pnt, if world2pix and pix2world worked correctly:'
print 'ra_out = ',ra_out,'  dec_out = ',dec_out  

# And now you want the value of the beam at these positions
selected_px, selected_py = w.wcs_world2pix(ra_withinradius, dec_withinradius, 1)

# Note for 300MHz radius size might cross over the beam edge, will need to filter out NaN values in this case.

# Filtering out sources outside of the beam image.
nan_ind = np.isnan(selected_py)==False

#print "length nan_ind = ",len(nan_ind),type(nan_ind),np.shape(nan_ind),np.shape(nan_ind)

selected_py = selected_py[np.isnan(selected_py)==False]
selected_px = selected_px[np.isnan(selected_px)==False]

print "Number of sources after NaN filtering = ",len(selected_py)

################################################################################
#Nan Filtering (Sources outside the beam)
################################################################################
name_withinradius=name_withinradius[nan_ind]
ra_withinradius=ra_withinradius[nan_ind]
dec_withinradius=dec_withinradius[nan_ind]
ra_str_withinradius=ra_str_withinradius[nan_ind]
dec_str_withinradius=dec_str_withinradius[nan_ind]
S200_fit_withinradius=S200_fit_withinradius[nan_ind]
alpha_withinradius=alpha_withinradius[nan_ind]
q_curve_withinradius=q_curve_withinradius[nan_ind]
#int_flux_subband_withinradius=int_flux_subband_withinradius[nan_ind]
int_flux_wide_withinradius=int_flux_wide_withinradius[nan_ind]
a_wide_withinradius=a_wide_withinradius[nan_ind]
b_wide_withinradius=b_wide_withinradius[nan_ind]
pa_wide_withinradius=pa_wide_withinradius[nan_ind]
peak_flux_wide_withinradius=peak_flux_wide_withinradius[nan_ind]
S_centralfreq_withinradius=S_centralfreq_withinradius[nan_ind]

# Converted to integers for indexing purposes.
pixel_x = np.rint(selected_px).astype('i')
pixel_y = np.rint(selected_py).astype('i')

# Vectorised beamvalue determination, no need for a for loop.
beamvalue = imagedata[0,0,pixel_y-1,pixel_x-1]

# Calculating the apparent flux density for each source within radius:
S_centralfreq_uncorrected = S_centralfreq_withinradius*beamvalue

# 0.0127 Jy is the theoretical 300 MHz sensitivity.
threshold = 0.0127
print 'threshold = ',threshold, ' Jy'

# Further filtering out nan values, sources outside of beam.
#
# Note: This filtering may be causing issues, look into later.

nan_ind2 = np.isnan(S_centralfreq_uncorrected)==False
S_centralfreq_uncorrected = S_centralfreq_uncorrected[nan_ind2]

print "Number of sources after secondary NaN filtering = ",len(S_centralfreq_uncorrected)

################################################################################
#Further nan and threshold filtering
################################################################################
keep_name=name_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_ra=ra_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_dec=dec_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_beamvalue=beamvalue[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_alpha=alpha_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_q=q_curve_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_S_centralfreq=S_centralfreq_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_S_centralfreq_uncorrected=S_centralfreq_uncorrected[S_centralfreq_uncorrected > threshold]
keep_ra_str=ra_str_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_dec_str=dec_str_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_int_flux_wide=int_flux_wide_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_a_wide=a_wide_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_b_wide=b_wide_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_pa_wide=pa_wide_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]
keep_peak_flux_wide=peak_flux_wide_withinradius[nan_ind2][S_centralfreq_uncorrected > threshold]

# These also need ot be filtered:
pixel_x = pixel_x[nan_ind2][S_centralfreq_uncorrected > threshold]
pixel_y = pixel_y[nan_ind2][S_centralfreq_uncorrected > threshold]

print "Number of sources after thresholding = ",len(pixel_x)

################################################################################
#Determining the beam spectral index and curvature coefficeints
################################################################################

#Initialising beam spectral index and curvature vectors.
alpha_beam = np.empty(len(pixel_x))
q_beam = np.empty(len(pixel_x))

print "Determining the beam spectral index and curvature"

start = time.time()
for i in range(len(keep_S_centralfreq_uncorrected)):
	Slogbeam_temp = np.log10(np.array(Beam_Cube[:,pixel_y[i]-1,pixel_x[i]-1]))
	if i==10:
		print "Example beam coefficients at pix_x = {0}, pix_y = {1}:".format(pixel_x[i]-1,pixel_y[i]-1)
		print Slogbeam_temp
	# Dividing by freq to normalise with respect to the centre of the band.
        C = Fiterr(Sqd_Quad_Residuals,[1,1,1],Nu/freq,Slogbeam_temp,1)
	alpha_beam[i] = C[1]
	q_beam[i] = C[2]

end = time.time()
print "Beam log-quadratic coefficient fitting duration = ", np.round(end - start,3),"s\n"

# Apparent spectral index and curvature for each source within radius.
keep_alpha_uncorrected = keep_alpha + alpha_beam
keep_q_uncorrected = keep_q + q_beam

print 'Number of PUMA sources selected = ',len(keep_S_centralfreq_uncorrected)
print 'Write these to the output file'

filename_body='model'

newvot = Table( [keep_name, keep_ra_str, keep_dec_str, keep_ra, keep_dec, keep_S_centralfreq, keep_alpha, alpha_beam, keep_alpha_uncorrected,keep_q, q_beam, keep_q_uncorrected, keep_beamvalue, keep_S_centralfreq_uncorrected, keep_peak_flux_wide, keep_int_flux_wide, keep_a_wide, keep_b_wide, keep_pa_wide], names=('Name','ra_str','dec_str','ra','dec','S_centralfreq','alpha',\
'alpha_beam','alpha_uncorrected','q','q_beam','q_uncorrected',\
'beamvalue','S_centralfreq_uncorrected','peak_flux_wide','int_flux_wide',\
'a_wide', 'b_wide', 'pa_wide'), meta={'name': 'first table'} )

writetoVO(newvot, filename_body+"_morecolumns_temp.vot")

# take only brightest 1500 sources

string1='stilts tpipe in='+filename_body+'_morecolumns_temp.vot cmd=\'sorthead -down 1500 S_centralfreq_uncorrected \' out='+filename_body+'_morecolumns.vot'

print string1

os.system(string1)

string2='stilts tpipe in='+filename_body+'_morecolumns_temp.vot omode=count | awk {{print $4}}'

string3='stilts tpipe in='+filename_body+'_morecolumns.vot omode=count | awk {{print $4}}'


print string2
print 'number of sources selected: '
os.system(string2)
print string3
print 'number of sources kept:'
os.system(string3)

end0 = time.time()
print "Total runtime = ", np.round(end0 - start0,3),"s\n"


