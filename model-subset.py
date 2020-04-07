#!/usr/bin/python

__author__ = "Jaiden Cook"
__credits__ = ["Jaiden Cook"]
__version__ = "1.0.1"
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

# Parser options:
from optparse import OptionParser

# Multiprocessing stuff:
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

# Scipy stuff:
import scipy
from scipy import stats

# casa-core stuff:
from casacore.tables import table,tablecolumn

# Astropy stuff:
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.io import fits
from astropy.io import ascii
from astropy.io.votable import parse_single_table
from astropy.table import Table,Column,vstack
from astropy.io.votable import writeto as writetoVO

def Beam_Read(Input):
    """
    Function reads in a beam, it then calls the Beam_Norm function to
    normalise that beam to the maximum sensitivity pixel. It then flattens
    that beam into a 1D vector which it then returns as the function output.
    """
    Beam = fits.open(Input)
    Beam_data = Beam[0].data
    Beam_map = Beam_data[0,0,:,:]

    return Beam_map,Beam

def Plot_Sky_Map(Map,Map_name,Map_label,scalemin,scalemax,plotsave=True):
    plt.clf()
    plt.imshow(Map,cmap='jet',origin='lower')
    plt.colorbar(label = '{0}'.format(Map_label))
    plt.xlabel("Pixels")
    plt.ylabel("Pixels")
    plt.clim(scalemin,scalemax)
    plt.tight_layout()
    if plotsave == True:
        plt.savefig("{0}".format(Map_name))
    else:
        plt.show()

if __name__ == "__main__":

	# Setting up parser options:
	usage="Usage: %prog [options]\n"
	parser = OptionParser(usage=usage)
	parser.add_option('--obsid',dest="obsid",default=None,help="Input OBSID for the observation, this should match the measurement set.")
	parser.add_option('--lobe',dest="lobe",default="PB",help="Input lobe ID (PB,GL1,GL2,GL3,GL4), the number of grating lobes (GL) can vary up to 4. Read the OBSID_lobe-metafile for lobe information.")
	parser.add_option('--scale',dest="scale",default=33.0,help="Pixel scale for imaging region, default is 33 asec which 2.5 times less than the 300 MHz predicted resolution.")
	parser.add_option('--subcond',dest="subcond",default=True,help="Boolean condition, defines whether the script subsets out the image region of interest, or writes image region model to measurement set.")
	parser.add_option('--v',dest="verbcond",default=False,help="Verbose output condition, default is false, useful for diagnostics.")

	(options, args) = parser.parse_args()

	obsid = int(options.obsid)
	condition = True

	# Setting the conditions.
	if str(options.subcond) == 'True':
		subcond = True
	elif str(options.subcond) == 'False':
		subcond = False
	else:
		print "subcond option needs to be set as True or False."

	if subcond == True:
		# This section works in conjunction with the python script lobe-finder.py, which outputs a dictionary containing the metadata for each
		# of the lobe regions. This part of the script takes a user inputted lobe region, generally predetermined by looking at the metadata
		# and determines the imaging parameters for WSCLEAN. This part of the script also subset the region of interest from the model image,
		# this then allows WSCLEAN to predict the visibilities of the rest of the sky, so the subsetted region can be imaged without contamination.

		if options.obsid == None:
			
			print "No --obsid given, see --h for help,terminating script."
			sys.exit(0)
		else:

			filename="{0}_deeper-model.fits".format(options.obsid)
			meta_file = open("{0}_lobe-metadata.txt".format(options.obsid),"r")

			if meta_file.mode == 'r':
				meta_dict = eval(meta_file.read())

			lobes = meta_dict['LOBE']

			index = 0
			for lobe in lobes:
				if lobe != str(options.lobe):
					index += 1
				elif lobe == str(options.lobe):
					break
				else:
					print "Lobe could not be found or does not match that stored in dictionary!"
					print "Quitting program"
					sys.exit(0)

			#"""
			Beammap,Beam = Beam_Read(filename)
			Beam_header = Beam[0].header
			
			# Get WCS for beam image
			w = wcs.WCS(Beam_header,naxis=2)
			
			# Getting the centre ra and dec for the given lobe.
			ra_cent = meta_dict['CENTRE-RA'][index]
			dec_cent = meta_dict['CENTRE-DEC'][index]
			
			# This is the size of the image.
			size = meta_dict['SIZE'][index]

			print 'ra_cent = {0} [deg]'.format(np.round(ra_cent,3)),'dec_cent = {0} [deg]'.format(np.round(dec_cent,3))
			print 'size = {0} [deg]'.format(np.round(2*size,3))

			pix_RA_cent, pix_DEC_cent = w.wcs_world2pix(ra_cent, dec_cent, 1)
			
			print pix_RA_cent,pix_DEC_cent

			ra_min = ra_cent - size
			ra_max = ra_cent + size
			dec_min = dec_cent - size
			dec_max = dec_cent + size

			# Wrapping the ra values back around.
			if ra_min < 0.0:
				ra_min = ra_min + 360.0
			if ra_max > 360.0:
				ra_max = ra_max - 360.0
			
			print "ra_min = {0} [deg]".format(np.round(ra_min,3)),"ra_max = {0} [deg]".format(np.round(ra_max,3))
			print "dec_min = {0} [deg]".format(np.round(dec_min,3)),"dec_max = {0} [deg]".format(np.round(dec_max,3))
			print "Delta RA = {0} [deg]".format(np.round(np.abs(ra_max-ra_min),3)), "Delta DEC = {0} [deg]".format(np.round(np.abs(dec_max-dec_min),3))

			# Determining the max and min x and y pixels.
			max_px, max_py = w.wcs_world2pix(ra_max, dec_max, 1)
			min_px, min_py = w.wcs_world2pix(ra_min, dec_min, 1)
			
			# Setting to integer scale.
			max_px = np.rint(max_px).astype('i') 
			min_px = np.rint(min_px).astype('i') 
			max_py = np.rint(max_py).astype('i') 
			min_py = np.rint(min_py).astype('i')

			print np.abs(max_px-min_px), np.abs(max_py-min_py)
			
			del_pix_x = np.abs(max_px-min_px)
			del_pix_y = np.abs(max_py-min_py)

			arcsec_pix_x = (3600*2*size)/del_pix_x
			arcsec_pix_y = (3600*2*size)/del_pix_y

			print "pix_x scale = {0} [asec/pix]".format(arcsec_pix_x)
			print "pix_y scale = {0} [asec/pix]".format(arcsec_pix_y)
			print "scale = {0}".format(options.scale)
			print "scale modifies = {0}".format(arcsec_pix_x/float(options.scale))

			# This is the scaled size of the new region image.
			img_x_pix_size = int(np.abs(max_px-min_px)*arcsec_pix_x/float(options.scale))
			img_y_pix_size = int(np.abs(max_py-min_py)*arcsec_pix_y/float(options.scale))

			print "Image centre: {0}".format(meta_dict['CENTRE-hms-dms'][index])
			print "New image dimensions {0}x{1}".format(img_x_pix_size, img_y_pix_size)

			#resolution = np.round(((2.5*options.scale)/60.0),1)
			resolution = np.round(((4.5*float(options.scale))/60.0),1)
			
			# Writing chgcentre and wsclean parameters to file.
			newfile = open("wsclean-{0}-{1}-par.txt".format(options.obsid,meta_dict["LOBE"][index]),"w")
			newfile.write("{0} {1} {2} {3} {4}".format(meta_dict['CENTRE-hms-dms'][index],img_x_pix_size,img_y_pix_size,options.scale,resolution))
			newfile.close()

			# This makes things easier.
			pix_range_x = [min_px,max_px]
			pix_range_y = [min_py,max_py]

			# Defining a new beam map.
			Beammap_nu = Beammap
			# Setting imaging region pixels to zero.
			#Beammap_nu[max_px:min_px,min_py:max_py] = 0.0
			# The x and the y axes are flipped in the beam dimensions, had issues with the incorrect regions being masked.
			Beammap_nu[min(pix_range_y):max(pix_range_y),min(pix_range_x):max(pix_range_x)] = 0.0
		
			# Defining the new beam map.
			nu_map = Beam
			nu_map[0].data[0,0,:,:] = Beammap_nu
			# Writing new image to the fits file.
			nu_map.writeto(filename,overwrite=True)
			#nu_map.writeto("Test.fits",overwrite=True)
			#"""
	elif subcond == False:

		# This options is after wsclean has predicted what the visibilities are for the rest of the sky (not including the subset specified in
		# in the previous section). This then gets subtracted from the corrected data and we should be left with the model of only the subsetted 
		# region.

		# This section is currently obselete.

		print "Separating out the subsetted data."

		t = table('{0}.ms'.format(obsid),readonly=False)
		
		# Loading in columns
		corrected_column = tablecolumn(t,'CORRECTED_DATA')
		model_column = tablecolumn(t,'MODEL_DATA')
	
		print np.mean(model_column[0])
		# Subtracting the sky model from the corrected data.
		model_column[0] = corrected_column[0] - model_column[0]
		print np.mean(model_column[0])
	
		# Saving data to measurement set.
		t.close()





