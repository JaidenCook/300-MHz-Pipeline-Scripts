# 300MHz MWA Observation Pipeline
This pipeline was developed during the duration of my Masters project to process and image MWA observations at 300MHz. The aim of this pipeline is to perform sky-model calibration on an MWA 300MHz observation using a sky-model at 300MHz (Total-300MHz-skymodel.tar.xz). This approach attenuates a model of the MWA tile beam to a subset of the sky-model to create an apparent sky-model. This is performed for a particular observation. This model is then used to calibrate the visibilities of an input observation using the software CALIBRATE (Offringa, 2016). Once the data is calibrated an all-sky image is produced using WSCLEAN (Offringa, 2014). This all-sky image is then used to self-calibrate the observation visibilities. Once the data has been further calibrated, a deep image of the primary beam is made by subtracting a model of the sky visiblities (minus the primary beam contribution) from the observation visibilites. This image is also produced with WSCLEAN. The software used in this pipeline is largely publically available, and is free.

## Pipeline Dependencies
The software used in this pipeline is listed below, and where possible a link is provided to download the required software.

### CASA
CALIBRATE (This is MWA software that is not necessarily publically available, https://github.com/ICRAR/mwa-reduce)

CASA-CORE python wrapper (http://casacore.github.io/casacore/ ,https://github.com/casacore/casacore)

CHGCENTRE (https://sourceforge.net/p/wsclean/wiki/Installation/)

COTTER (https://github.com/MWATelescope/cotter)

WSCLEAN (https://sourceforge.net/p/wsclean/wiki/Installation/)

### Python Packages

Aeagean (https://github.com/PaulHancock/Aegean/wiki/Installation-and-Requirements)

astropy

joblib

matplotlib

multiprocessing

mwa_client (https://github.com/ICRAR/manta-ray-client) - This is used to download MWA observations.

mwa_pb (https://github.com/MWATelescope/mwa_pb)

numpy

python version used is 2.7, may be updated to 3.6 at a later date.

python-casacore (https://github.com/casacore/python-casacore)

scipy

tqdm

**Note: A lot of these packages have their own dependencies, these are not listed.

# Pipeline Scripts

## Downloading Observations

Before observations are processed they need to be downloaded and pre-processed. This is done using the script 'obs-downloadV2.sh'. This script takes an input text file and csv file. The text file contains a list of observation ID's (OBSIDs) to be downloaded. The csv file contains a list of commands for each OBSID. These commands are passed to the mwa_client (refer to link in the previous section) which calls the ASVO server and submits the download jobs for each observation. This is all done within a python virtual environment. Once the raw files are downloaded they are consolidated into measurement sets by COTTER. The script flag-routine.sh is then called which applies RFI and baseline flagging methods from CASA. An example comman for this script is given below:

nohup obs-downloadV2.sh --input_dir=$MYDATA/Obs --obsid_list=$MYDATA/Obs/Obs_list.txt --csvfile=$MYDATA/Obs/Obs_dl_list.csv > obsdownload.log &

Obs_dl_list.csv is the input file required by the ASVO service to download observations and they should be of the format:

obs_id=1131038424, job_type=d, download_type=vis

It is reccomended to run this in the background since job times with the ASVO server are variable.

### Building the Apparent Sky-Model

Before the full pipeline script is discussed the individual scripts (which can be run independently) will be discussed. The first script builds the apparent sky-model which is used as a first pass calibration of the observation visibilities. This script is called build_appskyV4.sh. It takes only an input text file 'Obs_list.txt' which is the same one used for the observation download script. 

This script uses the 300 MHz sky-model. The path to this sky-model is hard coded into the script so this will need to be changed by the specific user. With this model the script Model_format_3.py to generate a model of the MWA tile beam for that observation. It generates this model for all the sources above the horizon in the sky-model. It then attenuates all these sources, and selects the brightest 1500 sources. It then creates chromatic model of the beam response across the 30.72MHz bandwidth for each of these sources. It outputs this as a vo.table which is then converted into a text file format (vo2newmodel_2.py) which can be read by CALIBRATE. An example of this input is given below:

nohup build_appskyV4.sh --input_dir=$MYDATA/Obs/${outputdir} --output_dir=$MYDATA/model/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > build-appsky-out.log &

# Calibrator Observations

This section is broken up into two categories. Calibrator observations and non-calibrator observations. The processing methods used for these two are slightly different.

### Sky Model Calibration

Once the observation apparent sky-model is produced it can be used to calibrate the observations visibilities. This is done by the script cal_year1_S300V2.sh. This has the same input text file. This script copies over the sky-model text file and runs CALIBRATE. Once the solutions are determined they are applied to the observation measurement set. Below is an example of how to run this script:

nohup cal_year1_S300V2.sh --input_dir=$MYDATA/Obs/${outputdir} --input_model=$MYDATA/model/${outputdir} --output_dir=$MYDATA/cal/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > cal-out.log &

### Self-Calibration

With the calibrated observation, granted the calibration solutions are good enough to image, self-calibration can occur. This is done with the script selfcal.sh. This script has the same input text file as the other scripts. This then calls WSCLEAN to make an all-sky image. The parameters are currently hard coded for this, but future versions may accept different parameters as inputs to the script. Once the all-sky image is made CALIBRATE is called again to perform the self-calibration.

nohup selfcal.sh --input_dir=$MYDATA/cal/${outputdir} --output_dir=$MYDATA/selfcal/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=$chan > selfcal-out.log &

# Non-calibrator Observations

These are observations that are not dominated by bright calibrator sources. In these situations the errors in the sky-model and the beam may lead to untenable calibration solutions. When this is the case calibration solutions may be transferred from an observation with valid calibration solutions. Granted this observation is at the same pointing and taken on the same night as the target observations. The only difference in the processing for these observations is that the calibration solutions from another observation need to be copied and tranferred to this one.

### Calibration Solution Transfer

The calibration transfer is performed by the script cal_transfer.sh. This script takes in the same input text file and the OBSID for the calibrator observation. This calibrator observation needs to have already been calibrated. It then copies over the solutions.bin file from the calibrator OBSID cal directory. These are then applied to the non-calibrator observations with the software applysolutions.

nohup cal_transfer.sh --input_dir=$MYDATA/Obs/${outputdir} --cal_dir=$MYDATA/cal/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --cal_obsid=${cal_obs} > cal_trans-out.log &

# Total Data Processing Pipeline

The total process described above is encapsulated into a single script S300-processing-pipeline.sh. This script does not contain the observation download and imagin scripts just the processing scripts. This takes in the same observation text file, but with a specific format. 

# Examples of Pipeline Inputs

## Obs_list.txt:

This is a text file which contains the list of obsid's to be processed, along with the tiles that need to be flagged and the observation directory. The format is specified as such:

OBSID TILE_LIST GROUP_ID CAL_OBSID

OBSID: The GPS time ID for the given observation.

TILE_LIST: sequence of integers that can range from 0-127, these are the ID's of the given tiles. If more than one tile is given then the ID's for each tile should be separated by a comma, i.e 0,105,80. If there are no tiles that need to be flagged then the string 'none' should be placed in the column.

GROUP_ID: This is the folder name for the observation group. Observations at 300 MHz using this strategy should come in groups of at least two OBSID's where the first observation is the calibrator observation. The calibrator observation is usually of a bright source such as Pictor A, Hydra A or 3C444. Hence the folder name is usually HydA{0}, where {0} is some integer that specifies the grid pointing of the MWA for that particular observation group.

CAL_OBSID: This is the OBSID of the calibrator observation, if the observation is a calibrator then this should just be left blank.

### Example of an Obs_list.txt File

1139931552 none HydA9 
1139935872 none HydA9 1139931552

or

1139931552 80,105 HydA9 
1139935872 80,105 HydA9 1139931552

**Note: This file only contains one observing group, but you could have multiple observing groups, this Obs_list.txt file should also be useable for the obs-download.sh script.

### Example Pipeline Command

nohup S300-processing-pipeline.sh --input_dir=$MYDATA --obsid_list=$MYDATA/Obs_list.txt --chan=236 > S300-pipeline-out.log &

S300-processing-pipeline.sh --input_dir=$MYDATA --obsid_list=$MYDATA/Obs_list.txt --chan=236

# Deep Imaging Script

The deep imaging is performed by the script deep-image.sh. This script performs source finding using Aegean on the all-sky image. It uses the output source catalogue to identify the grating sidelobes and the primary beam using the python script lobe-finder.py. This outputs a metadata text file in the format of a python dictionary. The meta data contains the centre, the size, the ID (PB,GL0,GL1...), and the number of sources per lobe. The input --lobeID can be used to select which lobe to image. Based on the input the model all-sky image output by WSCLEAN is then masked for that region using the information from the metafits file. This model image is then converted into a set of visibilities by WSCLEAN's predict function. These are then subtracted from the corrected visibilities. This process can be repeated N times. It is recommended that N <= 2. Each iteration requires all-sky imaging and model predicting. It becomes time consuming.

Once the subtraction is finished a deep image is created. The imaging parameters are hard coded. Future versions will allow for more user control over the imaging. For now these need to be manually changed. The deep image of the chosen lobe is then output. An example of the command is shown below:

deep-image.sh --input_dir=$MYDATA/selfcal/${outputdir} --output_dir=$MYDATA/image/${outputdir} --obsid_list=$MYDATA/Obs_list.txt --chan=236 --pipecond="no" --lobeID="PB"
