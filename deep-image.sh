#!/bin/bash

start_time=`date +%s`

# Set default values for optional parameters
mwapath="/home/jaidencook/.local/lib/python2.7/site-packages/mwa_pb/data/"
threshold=0.95
flag=yes
pipecond=no
lobeID=PB
Niter=1
flag_chans=4

# Read the options
TEMP=`getopt -o a:b:c:d:e:f:g:h: --long input_dir:,output_dir:,obsid_list:,chan:,pipecond:,lobeID:,Niter:,flag_chans: -- "$@"`
eval set -- "$TEMP"

# Extract options and their arguments into variables
while true ; do
      case "$1" in
        -a|--input_dir) # input directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) input_dir=$2 ; shift 2 ;;
            esac ;;
        -b|--output_dir) # output directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) output_dir=$2 ; shift 2 ;;
            esac ;;
        -c|--obsid_list) # obsID list (required argument)
            case "$2" in
                "") shift 2 ;;
                *) obsid_list=$2 ; shift 2 ;;
            esac ;;
        -d|--chan) # channel (required argument); set to 69, 93, 121, 145 or 169
            case "$2" in
                "") shift 2 ;;
                *) chan=$2 ; shift 2 ;;
            esac ;;
        -e|--pipecond) # channel (required argument); set to 69, 93, 121, 145 or 169
            case "$2" in
                "") shift 2 ;;
                *) pipecond=$2 ; shift 2 ;;
            esac ;;
        -f|--lobeID) # channel (required argument); set to 69, 93, 121, 145 or 169
            case "$2" in
                "") shift 2 ;;
                *) lobeID=$2 ; shift 2 ;;
            esac ;;
        -g|--Niter) # channel (required argument); set to 69, 93, 121, 145 or 169                                                      
            case "$2" in
                "") shift 2 ;;
                *) Niter=$2 ; shift 2 ;;
            esac ;;
        -h|--flag_chans) # input directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) flag_chans=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# Check required arguments have been specified
if [ -z "$input_dir" ]; then
  echo "Error: input_dir not specified. Aborting."
  exit 1
elif [ -z "$output_dir" ]; then
  echo "Error: output_dir not specified. Aborting."
  exit 1
elif [ -z "$obsid_list" ]; then
  echo "Error: obsid_list not specified. Aborting."
  exit 1
elif [ -z "$chan" ]; then
  echo "Error: chan not specified. Aborting."
  exit 1
fi

# Check chan parameter
if [ $chan != "69" ] && [ $chan != "93" ] && [ $chan != "121" ] && [ $chan != "145" ] && [ $chan != "169" ] && [ $chan != "236" ] && [ $chan != "237" ]; then
  echo "Error: chan should be set to 69, 93, 121, 145 or 169. Aborting."
  exit 1
fi

# Getting the OBSID.
if [ $pipecond == "yes" ]
then
  # In this case the input is a string not a file.
  obsid=$(echo $obsid_list | awk '{print $1}')
elif [ $pipecond == "no" ]
then
  # In this case the input is a file not a string.
  obsid=`awk '{print $1}' $obsid_list`
else
  echo "pipe condition not specified, exiting program."
  exit 1
fi

# Create output directory
if [ ! -e $output_dir ]; then
  mkdir $output_dir
fi

# Create snapshot directory
if [ -e $output_dir/$obsid ]; then
  rm -rf $output_dir/$obsid
fi
mkdir $output_dir/$obsid
cd $output_dir/$obsid

echo "======================================================================================="
echo "Coping over ${obsid}.* from directory ${input_data}/${obsid}"
echo "======================================================================================="

# Copy measurement set, metafits file and sky model to output directory
if [ -e $input_dir/$obsid/$obsid.ms ] && [ -e $input_dir/$obsid/$obsid.metafits ]; then
    cd $input_dir/$obsid/
    cd $output_dir/$obsid
    cp -r $input_dir/$obsid/$obsid.ms $input_dir/$obsid/${obsid}_lobe-metadata.txt $input_dir/$obsid/$obsid.ms.flagversions $input_dir/$obsid/$obsid.metafits  .
    cp -r $input_dir/$obsid/${obsid}_lobe-metadata.txt $input_dir/$obsid/${obsid}_deeper-model.fits .
else
    echo "Error: input files are missing. Aborting."
    exit 1
fi

echo "======================================================================================="
echo "Subsetting lobe : ${lobeID} from the rest of the sky."
echo "======================================================================================="

model-subset.py --obsid ${obsid} --lobe ${lobeID} --scale=18.0

# Scale of 18.0 is roughly 4.5 pixels per psf.

# Reading in wsclean deep image paramaters.
ra_hms=$(awk '{print $1}' wsclean-${obsid}-${lobeID}-par.txt)
dec_dms=$(awk '{print $2}' wsclean-${obsid}-${lobeID}-par.txt)
ra_img_size=$(awk '{print $3}' wsclean-${obsid}-${lobeID}-par.txt)
dec_img_size=$(awk '{print $4}' wsclean-${obsid}-${lobeID}-par.txt)
scale=$(awk '{print $5}' wsclean-${obsid}-${lobeID}-par.txt)
resolution=$(awk '{print $6}' wsclean-${obsid}-${lobeID}-par.txt)

echo "${lobeID} centre ${ra_hms}${dec_dms}"
echo "Image dimensions = ${ra_img_size}x${dec_img_size}"
echo "Pixel scale = ${scale} [asec]"
echo "resolution = ${resolution} [amin]"

echo "======================================================================================="
echo "Predicting residual sky visibilities."
echo "======================================================================================="

flag-routine.sh --obsid=${obsid} --data_column=corrected --flag_chans=${flag_chans}

# This is used for deeper imaging.
wsclean -predict -name ${obsid}_deeper -pol I -abs-mem 30 -j 12 $obsid.ms

echo "======================================================================================="
echo "Subtracting residual sky visibilities from the corrected data."
echo "======================================================================================="

# Suggested by Andre to updated corrected data column.
taql update ${obsid}.ms set CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA

#model-subset.py --obsid ${obsid} --subcond False

echo "======================================================================================="
echo "Performing $Niter round(s) of all-sky imaging and subtraction."
echo "======================================================================================="

# "Super Major" cleaning iterations.
Nseq=$(seq $Niter)

for i in $Nseq;
do
    wsclean -name ${obsid}_deeper -size 7000 7000 -niter 300000 -auto-threshold 1.0 -auto-mask 3.0 -pol I -weight uniform -scale 59asec -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.85 -minuv-l 60 -taper-gaussian 2.4amin $obsid.ms

    echo "======================================================================================="
    echo "Subsetting lobe : ${lobeID} from the rest of the sky."
    echo "======================================================================================="

    model-subset.py --obsid ${obsid} --lobe ${lobeID} --scale=18.0

    echo "======================================================================================="
    echo "Predicting residual sky visibilities."
    echo "======================================================================================="

    wsclean -predict -name ${obsid}_deeper -pol I -abs-mem 30 -j 12 $obsid.ms

    echo "======================================================================================="
    echo "Subtracting residual sky visibilities from the corrected data."
    echo "======================================================================================="

    taql update ${obsid}.ms set CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA

    #rm *.fits
    flag-routine.sh --obsid=${obsid} --data_column=corrected

done

##
#wsclean -name ${obsid}_deeper -size 7000 7000 -niter 300000 -auto-threshold 1.0 -auto-mask 3.0 -pol I -weight uniform \
#-scale 59asec -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.85 -minuv-l 60 -taper-gaussian 2.4amin $obsid.ms
#
#model-subset.py --obsid ${obsid} --lobe ${lobeID}
#
#wsclean -predict -name ${obsid}_deeper -pol I -abs-mem 30 -j 12 $obsid.ms
#
#taql update ${obsid}.ms set CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA
###
#
#rm *.fits

#flag-routine.sh --obsid=${obsid} --data_column="CORRECTED"

# The projection should only be changed after the prediction, since the image projection will be different to the data in the 
# measurement set.
chgcentre ${obsid}.ms ${ra_hms} ${dec_dms}

# It might be a good idea to make an all-sky image before making a deep image to see how well the subtraction works.

echo "======================================================================================="
echo "Imaging lobe : ${lobeID}"
echo "======================================================================================="

#wsclean -name ${obsid}_deeper -size 7000 7000 -niter 300000 -auto-threshold 1.0 -auto-mask 3.0 -pol I -weight uniform \
#-scale 59asec -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.85 -minuv-l 60 -taper-gaussian 2.4amin $obsid.ms

#wsclean -name ${obsid}_deeper -size ${ra_img_size} ${dec_img_size} -niter 300000 -auto-threshold 5.0 -auto-mask 8.0 -pol I -weight uniform \
#-scale "${scale}asec" -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.85 -minuv-l 60 -taper-gaussian "${resolution}amin" $obsid.ms

#wsclean -name ${obsid}_deeper -size ${ra_img_size} ${dec_img_size} -local-rms -niter 150 -nmiter 1 -auto-threshold 3.0 -auto-mask 5.0 -pol I -weight uniform \
#-scale "${scale}asec" -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.09 -minuv-l 60 -taper-gaussian "${resolution}amin" $obsid.ms

#wsclean -name ${obsid}_deeper -size ${ra_img_size} ${dec_img_size} -local-rms -niter 10000 -auto-threshold 1.0 -auto-mask 3.0 -pol I -weight uniform \
#-scale "${scale}asec" -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.5 -minuv-l 60 -taper-gaussian "${resolution}amin" $obsid.ms

# Silly me did this without using multi-scale, see how this changes the output.
wsclean -name ${obsid}_deep -multiscale -size ${ra_img_size} ${dec_img_size} -local-rms -niter 500000 -auto-threshold 1.0 -auto-mask 3.0 -pol I -weight briggs 0.0 \
-scale "${scale}asec" -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.5 -minuv-l 60  $obsid.ms

## Silly me did this without using multi-scale, see how this changes the output.
#wsclean -name ${obsid}_deep -multiscale -size ${ra_img_size} ${dec_img_size} -local-rms -niter 500000 -auto-threshold 1.0 -auto-mask 3.0 -pol I -weight uniform\
#-scale "${scale}asec" -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.5 -minuv-l 60 $obsid.ms

