#!/bin/bash

# Calibrate data using GLEAM year 1 catalogue
start_time=`date +%s`

# Set default values for optional parameters
flag=yes
pipecond=no

# Read the options
TEMP=`getopt -o a:b:c:d:e:f: --long input_dir:,output_dir:,cal_dir:,obsid_list:,cal_obsid:,pipecond: -- "$@"`
eval set -- "$TEMP"

# Extract options and their arguments into variables
while true ; do
      case "$1" in
        -a|--input_dir) # input directory containing raw measurement set (required argument)
            case "$2" in
                "") shift 2 ;;
                *) input_dir=$2 ; shift 2 ;;
            esac ;;
        -b|--output_dir) # output directory (required argument)                                         
            case "$2" in
                "") shift 2 ;;
                *) output_dir=$2 ; shift 2 ;;
            esac ;;
        -c|--cal_dir) # input directory containing calibrator solutions.
            case "$2" in
                "") shift 2 ;;
                *) cal_dir=$2 ; shift 2 ;;
            esac ;;
        -d|--obsid_list) # obsID list (required argument)
            case "$2" in
                "") shift 2 ;;
                *) obsid_list=$2 ; shift 2 ;;
            esac ;;
        -e|--cal_obsid) # obsID list (required argument)
            case "$2" in
                "") shift 2 ;;
                *) cal_obsid=$2 ; shift 2 ;;
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
elif [ -z "$cal_dir" ]; then
  echo "Error: cal_dir not specified. Aborting."
  exit 1
elif [ -z "$cal_dir" ]; then
  echo "Error: cal_dir not specified. Aborting."
  exit 1
elif [ -z "$obsid_list" ]; then
  echo "Error: obsid_list not specified. Aborting."
  exit 1
elif [ -z "$cal_obsid" ]; then
  echo "Error: No calibrator OBSID given, aborting script!"
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

# Group output directory should already exist if transferring calibration solutions.
if [ ! -e $cal_dir/$cal_obsid ]; then
   echo "Calibration obsid directory doesn't exits!"
   echo "Exiting program!"
   exit 0
fi
   
if [ ! -e $output_dir/$obsid ]; then
    mkdir $output_dir/$obsid
fi

cd $output_dir/$obsid

# Copy measurement set and the metafits file.
if [ -e $input_dir/$obsid/$obsid.ms ] && [ -e $input_dir/$obsid/$obsid.metafits ]; then
    cp -r $input_dir/$obsid/$obsid.ms $input_dir/$obsid/$obsid.metafits .
else
    echo "Error: input files are missing. Aborting."
    exit 1
fi

# -------------------------------------------------------------
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

################################################################################
#Calibration testing begin.
################################################################################

# Copying over the solutions bin from the calibrator observation. 
cp $cal_dir/${cal_obsid}/${cal_obsid}_solutions.bin .

# Apply the solutions
echo "Transferring calibration solutions from $cal_obsid to $obsid"
applysolutions $obsid.ms ${cal_obsid}_solutions.bin

# Flagging RFI.
flag-routine.sh --obsid=${obsid} --data_column=corrected

## Further flagging of the corrected data is required to eliminate RFI.                                                                  
#echo "======================================================================================="
#echo "Baseline and RFI flagging."
#echo "======================================================================================="
#
## Performing further RFI flagging.                                                                                                      
#echo "Applying rflag."
#echo "======================================================================================="
#
#casa --nologger -c "flagdata('${obsid}.ms',mode='rflag',datacolumn='CORRECTED',timedevscale=4.0,action='apply')"
#
#echo "======================================================================================="
#echo "Applying tcrop."
#echo "======================================================================================="
#
#casa --nologger -c "flagdata('${obsid}.ms',mode='tfcrop',datacolumn='CORRECTED',freqcutoff=3.0,usewindowstats='sum',action='apply')"
#
#echo "======================================================================================="
#echo "Antennas to be flagged."
#echo "======================================================================================="
#
## Getting a list of antenna baselines.                                                                                                  
#ANT_FLAG=$(ms_flag_by_uvdist.py "${obsid}.ms" "CORRECTED_DATA")
#echo "Antenna pairs ot be flagged."
#echo $ANT_FLAG
#
#echo "======================================================================================="
#echo "Applying steflag."
#echo "======================================================================================="
#
## Flagging baselines.                                                                                                                   
#casa --nologger -c "flagdata(vis='${obsid}.ms',antenna='$ANT_FLAG',mode='manual',action='apply')"

echo "======================================================================================="
echo "Plotting CORRECTED_DATA channels and baselines."
echo "======================================================================================="

casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='frequency',yaxis='amp',correlation='xx,yy',ydatacolumn='corrected',coloraxis='spw',plotfile='amp-v-freq_{0}_corrected.png'.format(${obsid}),showgui=False,overwrite=True)"
casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='baseline',yaxis='amp',correlation='xx,yy',ydatacolumn='corrected',coloraxis='spw',plotfile='amp-v-base_{0}_corrected.png'.format(${obsid}),showgui=False,overwrite=True)"

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
