# 300-MHz-Pipeline-Scripts
This repository contains all the scripts from my honours degree that I used to calbrate 300 MHz observations. These scripts are the modified calibration and imaging pipeline described in (https://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/GLEAM_year2), where they have been augmented to operate at 300 MHz.

## Modifications to the Original Pipeline:
The scripts build_beamsV2.sh and build_appskyV2.sh replace the build_model.sh script from the original pipeline. In future these two scripts will be consolidated again, since the reason they were split into two is no longer an issue.

There are two versions of attenuate_by_beam.py. This script is where most of the modifications have occured, it has been modified from the original script to calculate the curvature in the beam and hence the apparent curvature for sources in the sky model. Additionally the order of operations has been changed and unnecessary for loops have been modified to make the script run more efficiently, this is necessary when fitting for beam curvature. There is a 300 MHz version which operates using the 300 MHz sky model that I developed for my honours thesis, if you are interested in this model email me (jaiden.cook@student.curtin.edu.au). The second version uses the full as yet releases GLEAM catalogue, and has been further generalised to determine the beam curvature for any specified MWA frequency band. The new GLEAM version of attenuate_by_beam_V7.py should work with the old build_model.sh, though the inputs parsed to the script have now changed, so build_model.sh will need a slight modifcation when calling attenuate_by_beam_V7.py, below is the suggested change:

$srun attenuate_by_beam_S300V7.py --pbeam ${obsid}-MFS-psf_beamI.fits --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold

Important note: the threshold in the current version of attenuate_by_beam_S300V7.py and the GLEAM variant is manually set in the code, this will need ot be commented out if you want to be able to parse your own threholds, or manually changed.

Alternatively the GLEAM variant should work with build_appskyV2.sh, simply the name will need to be changed. I have only tested that this works with observations at 300 MHz, so if there are any issues please send me an email.

The other scripts have simply been modifed to allow for the 300 MHz channel, apart from this change, they are effectively the same as the original scripts.
