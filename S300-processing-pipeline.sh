#!/bin/bash

# Calibrate data using GLEAM year 1 catalogue
start_time=`date +%s`

# Future potential options will include cpus, and mem.

# Read the options
TEMP=`getopt -o a:b:c: --long input_dir:,obsid_list:,chan: -- "$@"`
eval set -- "$TEMP"

# Extract options and their arguments into variables
while true ; do
      case "$1" in
        -a|--input_dir) # input directory containing model for calibration (required argument)
            case "$2" in
                "") shift 2 ;;
                *) input_dir=$2 ; shift 2 ;;
            esac ;;
        -b|--obsid_list) # obsID list (required argument)
            case "$2" in
                "") shift 2 ;;
                *) obsid_list=$2 ; shift 2 ;;
            esac ;;
        -c|--chan) # channel (required argument); set to 69, 93, 121, 145 or 169
            case "$2" in
                "") shift 2 ;;
                *) chan=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# The number of observations.                                                                        
N_lines=`wc -l $obsid_list | awk '{print $1}'`
echo "Number of observations = ${N_lines}"

# Setting the sequence for the for loop.                                                             
Nseq="$(seq -s ' ' 1 ${N_lines})"

for i in $Nseq
do
    line_number=$i
    # ith obsid.                                                                                     
    obslist_row=`awk "NR==${line_number}{print $1}" $obsid_list`
    obsid=`awk "NR==${line_number}{print $1}" $obsid_list | awk '{print $1}'`
    tiles=`awk "NR==${line_number}{print $1}" $obsid_list | awk '{print $2}'`
    outputdir=`awk "NR==${line_number}{print $1}" $obsid_list | awk '{print $3}'`
    cal_obs=`awk "NR==${line_number}{print $1}" $obsid_list | awk '{print $4}'`

    if [ -z $cal_obs ]
    then
        echo "======================================================================================="
        echo "#${line_number}"
        echo "Observation group: ${outputdir}"
        echo "OBSID: ${obsid}"
        echo "Flagged Tiles : ${tiles}"
        echo "======================================================================================="
    else
        echo "======================================================================================="
        echo "#${line_number}"
        echo "Observation group: ${outputdir}"
        echo "OBSID: ${obsid}"
        echo "Flagged Tiles : ${tiles}"
        echo "Calibrator OBSID : ${cal_obs}"
        echo "======================================================================================="
    fi

    echo "======================================================================================="
    echo "Building sky model for ${obsid}"
    echo "======================================================================================="


	# Building the apparent sky model.
    build_appskyV4.sh --input_dir=${input_dir}/Obs/${outputdir} --output_dir=${input_dir}/model/${outputdir} \
     --obsid_list="${obslist_row}" --chan=$chan --pipecond="yes"

    if [ -z $cal_obs ]
    then
    	# Processing case for calibrator observations.
    	# Initial calibration with the apparent sky model and further flagging of corrected data.
        echo "======================================================================================="
        echo "Calibrating apparent sky model for ${obsid}"
        echo "======================================================================================="
    	
        cal_year1_S300V2.sh --input_dir=${input_dir}/Obs/${outputdir} --input_model=${input_dir}/model/${outputdir}\
         --output_dir=${input_dir}/cal/${outputdir} --obsid_list="${obslist_row}" --chan=$chan --pipecond="yes"

        echo "======================================================================================="
        echo "Performing self-calibration on ${obsid}"
        echo "======================================================================================="

        # Self-calibration.
        selfcal.sh --input_dir=${input_dir}/cal/${outputdir} --output_dir=${input_dir}/selfcal/${outputdir}\
         --obsid_list="${obslist_row}" --chan=$chan --pipecond="yes"
    else
        echo "======================================================================================="
        echo "Transferring calibration solutions from ${cal_obs} to ${obsid}"
        echo "======================================================================================="

    	# Processing case for non-calibrator observations.
    	# Transferring calibration solutions from calibrator observation.
	    cal_transfer.sh --input_dir=${input_dir}/Obs/${outputdir} --output_dir=${input_dir}/cal/${outputdir}\
	    --cal_dir=${input_dir}/selfcal/${outputdir} --obsid_list="${obslist_row}" --cal_obsid=${cal_obs} --pipecond="yes"
    
        echo "======================================================================================="
        echo "Performing all-sky imaging ${obsid}"
        echo "======================================================================================="

        # Self-calibration.
        selfcal.sh --input_dir=${input_dir}/cal/${outputdir} --output_dir=${input_dir}/selfcal/${outputdir}\
         --obsid_list="${obslist_row}" --chan=$chan --pipecond="yes" --selfcal="no"

        # It is wise to be cautious when performing self-calibration on non-calibrator observations.

    fi
    

    # Self-calibration.
    #selfcal.sh --input_dir=${input_dir}/cal/${outputdir} --output_dir=${input_dir}/selfcal/${outputdir}\
    # --obsid_list="${obslist_row}" --chan=$chan --pipecond="yes" --selfcal="no"

done



