#!/bin/bash

start_time=`date +%s`

# Set default values for optional parameters
mwapath="/home/jaidencook/.local/lib/python2.7/site-packages/mwa_pb/data/"
threshold=0.95
flag=yes
pipecond=no
selfcal=yes

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
        -f|--selfcal) # channel (required argument); set to 69, 93, 121, 145 or 169
            case "$2" in
                "") shift 2 ;;
                *) selfcal=$2 ; shift 2 ;;
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

echo "========================================================================"
echo "Coping over ${obsid}.* from directory ${input_data}/${obsid}"
echo "========================================================================"

# Copy measurement set, metafits file and sky model to output directory
if [ -e $input_dir/$obsid/$obsid.ms ] && [ -e $input_dir/$obsid/$obsid.metafits ]; then
    cd $input_dir/$obsid/
    cd $output_dir/$obsid
    cp -r $input_dir/$obsid/$obsid.ms $input_dir/$obsid/$obsid.ms.flagversions $input_dir/$obsid/$obsid.metafits .
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

maxuvm=887250/$chan

# Changing the centre to zenith.
chgcentre -zenith ${obsid}.ms

echo "======================================================================================="
echo "Making a shallow all-sky image."
echo "======================================================================================="

# Creating a shallow all sky image with 4arc minute resolution. This is so we can identify bright regions.
# Parameters for the all-sky image are fixed for the case of phase 1 300MHz.
wsclean -name ${obsid}_deeper -size 7000 7000 -niter 300000 -auto-threshold 0.3 -auto-mask 3.0 -pol I -weight uniform \
-scale 59asec -abs-mem 30 -j 12 -apply-primary-beam -mwa-path $mwapath -mgain 0.85 -minuv-l 60 -taper-gaussian 2.4amin $obsid.ms

echo "======================================================================================="
echo "Running BANE and Aegean."
echo "======================================================================================="

# After the all-sky image is created Aegean is run to determine where the flux density is.
BANE ${obsid}_deeper-image.fits
aegean --autoload --table ${obsid}-sources.fits ${obsid}_deeper-image.fits

echo "======================================================================================="
echo "Identifying primary beam and grating lobes."
echo "======================================================================================="

lobe-finder.py --obsid ${obsid}

echo "======================================================================================="
echo "Performing self-calibration"
echo "======================================================================================="

# For calibration transfer observations we may not want to self-calibrate.
if [ $selfcal == "yes" ]
then
  # Performing self cal:
  calibrate -absmem 30 -minuv 60 $obsid.ms ${obsid}_solutions.bin
  
  echo "======================================================================================="
  echo "Applying solutions to $obsid"
  echo "======================================================================================="
  
  # Apply the solutions
  applysolutions $obsid.ms ${obsid}_solutions.bin
  
  # Plotting the amplitude of the channels vs frequency.
  echo "======================================================================================="
  echo "Creating phase and amplitude plots."
  echo "======================================================================================="
  
  # Plot phase and amplitude calibration solutions
  aocal_plot.py --refant=127 ${obsid}_solutions.bin
else
  echo "Not selfcalibrating"
fi

end_time=`date +%s`
duration=`echo "$end_time-$start_time" | bc -l`
echo "Total runtime = $duration sec"

exit 0


