
## Pipeline Dependencies: (I will need to give links as to where to get these)

# CASA
CASA-CORE python wrapper
CHGCENTRE
COTTER 
WSCLEAN - These need casa-core to function.

python packages:
python version used is 2.7, may be updated to 3.6 at a later date.
astropy
matplotlib
mwa_pb (github download for this)
numpy
scipy
joblib
multiprocessing
tqdm

Note: If processing on Pawsey scripts will need to be modified.
Self-Note: I will also have to give the version of all the packages I used.

## Downloading observations:

nohup obs-downloadV2.sh --input_dir=$MYDATA/Obs --obsid_list=$MYDATA/Obs/Obs_list.txt --csvfile=$MYDATA/Obs/Obs_dl_list.csv > obsdownload.log &

Obs_dl_list.csv is the input file required by the ASVO service to download observations and they should be of the format:

obs_id=1131038424, job_type=d, download_type=vis

COTTER is manually applied to the gpu.fits files of the observation as well as manual flagging using CASA.

## Total data processing pipeline:
"
This does not include the observation downloading and inital flagging or the image script.
"
nohup S300-processing-pipeline.sh --input_dir=$MYDATA --obsid_list=$MYDATA/Obs_list.txt --chan=236 > S300-pipeline-out.log &

S300-processing-pipeline.sh --input_dir=$MYDATA --obsid_list=$MYDATA/Obs_list.txt --chan=236

### Calibrator Observations:

## Individual scripts:
#
# Building the apparent sky model:
nohup build_appskyV4.sh --input_dir=$MYDATA/Obs/${outputdir} --output_dir=$MYDATA/model/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > build-appsky-out.log &

# Sky model calibration:
nohup cal_year1_S300V2.sh --input_dir=$MYDATA/Obs/${outputdir} --input_model=$MYDATA/model/${outputdir} --output_dir=$MYDATA/cal/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > cal-out.log &

# Self-calibration:
nohup selfcal.sh --input_dir=$MYDATA/cal/${outputdir} --output_dir=$MYDATA/selfcal/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > selfcal-out.log &

### Non-calibrator Observations:

## Individual scripts:
#
# Building the apparent sky model:(This is optional in this case)
nohup build_appskyV4.sh --input_dir=$MYDATA/Obs/${outputdir} --output_dir=$MYDATA/model/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > build-appsky-out.log &

# Calibration solution transfer:
nohup cal_transfer.sh --input_dir=$MYDATA/Obs/${outputdir} --cal_dir=$MYDATA/cal/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --cal_obsid=${cal_obs} > cal_trans-out.log &

# Self-calibration:
nohup selfcal.sh --input_dir=$MYDATA/cal/${outputdir} --output_dir=$MYDATA/selfcal/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > selfcal-out.log &

## Inputs:

# Obs_list.txt:

This is a text file which contains the list of obsid's to be processed, along with the tiles that need to be flagged and the observation directory. The format is specified as such:

OBSID TILE_LIST GROUP_ID CAL_OBSID

OBSID: The GPS time ID for the given observation.

TILE_LIST: sequence of integers that can range from 0-127, these are the ID's of the given tiles. If more than one tile is given then the ID's for each tile should be separated by a comma, i.e 0,105,80. If there are no tiles that need to be flagged then the string 'none' should be placed in the column.

GROUP_ID: This is the folder name for the observation group. Observations at 300 MHz using this strategy should come in groups of at least two OBSID's where the first observation is the calibrator observation. The calibrator observation is usually of a bright source such as Pictor A, Hydra A or 3C444. Hence the folder name is usually HydA{0}, where {0} is some integer that specifies the grid pointing of the MWA for that particular observation group.

CAL_OBSID: This is the OBSID of the calibrator observation, if the observation is a calibrator then this should just be left blank.

# Example of an Obs_list.txt file:

1139931552 none HydA9 
1139935872 none HydA9 1139931552

or

1139931552 80,105 HydA9 
1139935872 80,105 HydA9 1139931552

**Note: This file only contains one observing group, but you could have multiple observing groups, this Obs_list.txt file should also be useable for the obs-download.sh script.

For processing outside of the pipeline, the Obs_list.txt file requires only one row entry. This can be a calibrator or non-calibrator, just note that a non-calibrator will require a calibrated source for processing.

## Deep Imaging Script:

deep-image.sh --input_dir=$MYDATA/selfcal/${outputdir} --output_dir=$MYDATA/image/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=236 --pipecond="no" --lobeID="PB"
