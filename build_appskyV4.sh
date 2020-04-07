#!/bin/bash

# Build sky model for snapshot using the 300 MHz interpolated catalogue.

start_time=`date +%s`

# Set default values for optional parameters
resolution=0.022 # Approximate 300 MHz psf major and minor axis size.
pipecond="no"

# Read the options
TEMP=`getopt -o a:b:c:d:e: --long input_dir:,output_dir:,obsid_list:,chan:,pipecond: -- "$@"`
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
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# Check required arguments have been specified
if [ -z "$input_dir" ]; then
  echo "Error: input_dir not specified. Aborting."
  exit 1
elif [ -z "$obsid_list" ]; then
  echo "Error: obsid_list not specified. Aborting."
  exit 1
fi

# Check chan parameter
if [ $chan != "69" ] && [ $chan != "93" ] && [ $chan != "121" ] && [ $chan != "145" ] && [ $chan != "169" ] && [ $chan != "236" ] && [ $chan != "237" ]; then
  echo "Error: chan should be set to 69, 93, 121, 145, 169 or 236. Aborting."
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

# Create snapshot directory
if [ -e $output_dir ]; then
  echo "Error: Output directory $output_dir already exists. Aborting."
  #exit 1
else
  mkdir $output_dir
fi

# Create snapshot directory
if [ -e $output_dir/$obsid ]; then
  echo "Error: Output directory $output_dir/$obsid already exists. Aborting."
  exit 1
else
  mkdir $output_dir/${obsid}
  cd $output_dir/$obsid
fi

# Coping over the obsid metafits file.
if [ -e $input_dir/$obsid/$obsid.metafits ]; then 
  cp $input_dir/$obsid/$obsid.metafits .
else
  echo "Error: input files are missing. Aborting."
  exit 1
fi

# Set delays and frequency.
delays=$(fitshdr $obsid.metafits | grep 'DELAYS' | awk '{print $3}')
freq=$(fitshdr $obsid.metafits | grep 'FREQCENT' | awk '{print $2}')

echo "$obsid measurement set and metafits copied to output directory."
echo "delays = $delays"

# Set frequency depending on channel
if [[ $chan -eq 69 ]]; then
  freq=88
elif [[ $chan -eq 93 ]]; then
  freq=118
elif [[ $chan -eq 121 ]]; then
  freq=154
elif [[ $chan -eq 145 ]]; then
  freq=185
elif [[ $chan -eq 169 ]]; then
  freq=215
elif [[ $chan -eq 236 ]]; then
  freq=300
elif [[ $chan -eq 237 ]]; then
  freq=303
fi
echo "Central frequency for channel $chan = $freq MHz"

# Get RA and Dec of pointing centre from .metafits file, in deg
ra_pnt=$(fitshdr $obsid.metafits | grep 'RA of pointing center' | awk '{print $3}')
dec_pnt=$(fitshdr $obsid.metafits | grep 'Dec of pointing center' | awk '{print $3}')

# Sky model at 300 MHz
catalogue=/home/jaidencook/Documents/Masters/catalogues/Total-300MHz-skymodel.fits

# Creating the output V02 table skymodel.
#Model_format_2.py --obsid $obsid --freq $freq --delays $delays --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold
Model_format_3.py --obsid $obsid --freq $freq --delays $delays --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt

# Convert the edited catalogue to a format which is readable by calibrate.
# Sources with int_flux_wide/peak_flux_wide < resolution will be considered to be unresolved.
vo2newmodel_2.py --catalogue model_morecolumns_temp.vot --output skymodelformat.txt --freq $freq --fluxcol S_centralfreq_uncorrected --coeff apparent_poly_coeff \
--acol a_wide --bcol b_wide --pacol pa_wide --point --resolution=$resolution

# -------------------------------------------------------------

end_time=`date +%s`
duration=`echo "$end_time-$start_time" | bc -l`
echo "Total runtime = $duration sec"

exit 0

