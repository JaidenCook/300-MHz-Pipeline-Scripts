#!/bin/bash

# Calibrate data using GLEAM year 1 catalogue
start_time=`date +%s`

# Set default values for optional parameters
flag=yes
pipecond=no

# Read the options
TEMP=`getopt -o a:b:c:d:e:f: --long input_dir:,input_model:,output_dir:,obsid_list:,chan:,pipecond: -- "$@"`
eval set -- "$TEMP"

# Extract options and their arguments into variables
while true ; do
      case "$1" in
        -a|--input_dir) # input directory containing raw measurement set (required argument)
            case "$2" in
                "") shift 2 ;;
                *) input_dir=$2 ; shift 2 ;;
            esac ;;
        -b|--input_model) # input directory containing model for calibration (required argument)
            case "$2" in
                "") shift 2 ;;
                *) input_model=$2 ; shift 2 ;;
            esac ;;
        -c|--output_dir) # output directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) output_dir=$2 ; shift 2 ;;
            esac ;;
        -d|--obsid_list) # obsID list (required argument)
            case "$2" in
                "") shift 2 ;;
                *) obsid_list=$2 ; shift 2 ;;
            esac ;;
        -e|--chan) # channel (required argument); set to 69, 93, 121, 145 or 169
            case "$2" in
                "") shift 2 ;;
                *) chan=$2 ; shift 2 ;;
            esac ;;
        -f|--pipecond) # channel (required argument); set to 69, 93, 121, 145 or 169
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
elif [ -z "$input_model" ]; then
  echo "Error: input_model not specified. Aborting."
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
if [ $chan != "69" ] && [ $chan != "93" ] && [ $chan != "121" ] && [ $chan != "145" ] && [ $chan != "169" ] && [ $chan != "236" ]; then
  echo "Error: chan should be set to 69, 93, 121, 145, 169, or 236. Aborting."
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

# Set other input parameters
ncpus=12
mem=31

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
echo "Coping over ${obsid}.*"
echo "======================================================================================="

# Copy measurement set, metafits file and sky model to output directory
if [ -e $input_dir/$obsid/$obsid.ms ] && [ -e $input_dir/$obsid/$obsid.metafits ] && [ -e $input_model/$obsid/skymodelformat.txt ]; then
# kluge to refresh stale file handles
    cd $input_dir/$obsid/
    cd $output_dir/$obsid
    cp -r $input_dir/$obsid/$obsid.ms $input_dir/$obsid/$obsid.ms.flagversions $input_dir/$obsid/$obsid.metafits $input_model/$obsid/skymodelformat.txt .
else
    echo "Error: input files are missing. Aborting."
    exit 1
fi

if [ $pipecond == "yes" ]
then
  # In this case the input is a string not a file.
  # Flag tiles if required
  if [ $flag == "yes" ]; then
    tile_list=$(echo $obsid_list | awk '{print $2}')
    echo $tile_list
    if [ -z "$tile_list" ]; then
      echo "No tiles to flag for snapshot $obsid"
    elif [ $tile_list == "none" ]; then
      echo "No tiles to flag for snapshot $obsid"
    else
      tile_list=`echo ${tile_list//,/ }`
  # Check tiles are integers between 0 and 127
      for tile in ${tile_list[@]}; do  
        if [ "${tile//[0-9]}" != "" ] || [ $(echo "$tile < 0"|bc) -eq 1 ] || [ $(echo "$tile > 127"|bc) -eq 1 ]; then
          echo "Error: tile $tile is not an integer between 0 and 127. Aborting."
          exit 1
        fi
      done
  # Flag tiles
      echo "Flagging tiles $tile_list for snapshot $obsid listed in $obsid_list"
      flagantennae $obsid.ms $tile_list
    fi
  fi
elif [ $pipecond == "no" ]
then
  # In this case the input is a file not a string.
  if [ $flag == "yes" ]; then
    tile_list=` awk '{print $2}' $obsid_list`
    echo $tile_list
    if [ -z "$tile_list" ]; then
      echo "No tiles to flag for snapshot $obsid"
    elif [ $tile_list == "none" ]; then
      echo "No tiles to flag for snapshot $obsid"
    else
      tile_list=`echo ${tile_list//,/ }`
  # Check tiles are integers between 0 and 127
      for tile in ${tile_list[@]}; do  
        if [ "${tile//[0-9]}" != "" ] || [ $(echo "$tile < 0"|bc) -eq 1 ] || [ $(echo "$tile > 127"|bc) -eq 1 ]; then
          echo "Error: tile $tile is not an integer between 0 and 127. Aborting."
          exit 1
        fi
      done
  # Flag tiles
      echo "Flagging tiles $tile_list for snapshot $obsid listed in $obsid_list"
      flagantennae $obsid.ms $tile_list
    fi
  fi
else
  echo "pipe condition not specified, exiting program."
  exit 1
fi

# determin max uv
maxuvm=887250/$chan

# This method takes the apparent sky model and calibrates the data.
echo "======================================================================================="
echo "Running CALIBRATE on ${obsid}."
echo "======================================================================================="

calibrate -m skymodelformat.txt -absmem $mem -minuv 60 -maxuv $maxuvm $obsid.ms ${obsid}_solutions.bin

echo "======================================================================================="
echo "Applying calibration solutions."
echo "======================================================================================="

applysolutions $obsid.ms ${obsid}_solutions.bin

casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='frequency',yaxis='amp',correlation='xx,yy',ydatacolumn='corrected',coloraxis='spw',plotfile='amp-v-freq_${obsid}_corrected.png',showgui=False,overwrite=True)"
casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='baseline',yaxis='amp',correlation='xx,yy',ydatacolumn='corrected',coloraxis='spw',plotfile='amp-v-base_${obsid}_corrected.png',showgui=False,overwrite=True)"

# Flagging the corrected data for RFI.
flag-routine.sh --obsid=${obsid} --data_column=corrected

casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='frequency',yaxis='amp',correlation='xx,yy',ydatacolumn='corrected',coloraxis='spw',plotfile='amp-v-freq_${obsid}_corrected-flgd.png',showgui=False,overwrite=True)"
casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='baseline',yaxis='amp',correlation='xx,yy',ydatacolumn='corrected',coloraxis='spw',plotfile='amp-v-base_${obsid}_corrected-flgd.png',showgui=False,overwrite=True)"

# Plot phase and amplitude calibration solutions
echo "Doing plot of phase and amplitude"
aocal_plot.py --refant=127 ${obsid}_solutions.bin

# Re-plot amplitude calibration solutions, this time setting maximum of y axis to 100 and 10000
for amp in 100 10000; do
  mkdir t
  aocal_plot.py --refant=127 --outdir=./t --amp_max=$amp ${obsid}_solutions.bin
  mv ./t/${obsid}_solutions_amp.png ${obsid}_solutions_amp_max${amp}.png
  rm -rf t
done

end_time=`date +%s`
duration=`echo "$end_time-$start_time" | bc -l`
echo "Total runtime = $duration sec"

exit 0
