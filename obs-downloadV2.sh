#!/bin/bash
start_time=`date +%s`

# Read the options
TEMP=`getopt -o a:b:c: --long input_dir:,obsid_list:,csvfile: -- "$@"`
eval set -- "$TEMP"

# Extract options and their arguments into variables
while true ; do
      case "$1" in
        -a|--input_dir) # input directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) input_dir=$2 ; shift 2 ;;
            esac ;;
        -b|--obsid_list) # obsID list (required argument)
            case "$2" in
                "") shift 2 ;;
                *) obsid_list=$2 ; shift 2 ;;
            esac ;;
        -c|--csvfile) # obsID list (required argument)
            case "$2" in
                "") shift 2 ;;
                *) csvfile=$2 ; shift 2 ;;
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

#source "/home/jaidencook/env/bin/activate"

# Downloading the observations.
mwa_client --csv="${csvfile}" --dir=$input_dir

# Getting out ofthe virtual environment, this is necessary otherwise cotter/casa wont work.
"deactivate"

for i in $Nseq
do
    line_number=$i
    # ith obsid.
    obsid=`awk "NR==${line_number}{print $1}" $obsid_list | awk '{print $1}'`
    outputdir=`awk "NR==${line_number}{print $1}" $obsid_list | awk '{print $3}'`
    
    # This is for basic testing purposes.
    echo "========================================================================"
    echo "#${line_number}                                                       ||"
    echo "OBSID:${obsid}                                                        ||"
    echo "Output directory = ./${outputdir}                                     ||"
    echo "========================================================================"

    # If the output directory doesn't exist make it.
    if [ ! -d ./${outputdir} ]; then
        mkdir -p ./${outputdir};
    fi
    
    # Move all files to the output directory, move to that directory, and unzip the gpu.fits files.
    mv ${obsid}* ./${outputdir}
    cd ./${outputdir}

    # If the obsid directory doesn't exist make it. It shouldn't exist.
    if [ ! -d ./${obsid} ]; then
        mkdir -p ./${obsid};
    fi

    mv ${obsid}* ./${obsid}
    # Move to new directory named after the obsid.
    cd ./${obsid}
    
    # Unzip gpu.fits files for cottering.
    unzip ${obsid}_vis.zip

    # Applying cotter.
    cotter -j 6  -mem 30 -timeres 4 -freqres 40 -edgewidth 80 -allowmissing -flagsubband 0,1,2 -flagdcchannels -m ${obsid}.metafits -o ${obsid}.ms ${obsid}*gpubox*.fits

    flag-routine.sh --obsid=${obsid} --data_column=data

    #echo "========================================================================"
    #echo "|| Baseline and RFI flagging.                                         ||"
    #echo "========================================================================"
    #
    ## Performing further RFI flagging.
    #echo "|| Applying rflag.                                                    ||"
    #echo "========================================================================"
    #
    #casa --nologger -c "flagdata('${obsid}.ms',mode='rflag',datacolumn='DATA',timedevscale=4.0,action='apply')"
    #
    #echo "========================================================================"
    #echo "|| Applying tcrop.                                                    ||"
    #echo "========================================================================"
    #
    #casa --nologger -c "flagdata('${obsid}.ms',mode='tfcrop',datacolumn='DATA',freqcutoff=3.0,usewindowstats='sum',action='apply')"

    ## Getting a list of antenna baselines.
    #ANT_FLAG=$(ms_flag_by_uvdist.py "${obsid}.ms" "DATA")
    #
    #echo "========================================================================"
    #echo "|| Antennas to be flagged.                                            ||"
    #echo "========================================================================"
    
    #echo $ANT_FLAG

    ## Flagging baselines.
    #echo "========================================================================"
    #echo "|| Applying steflag.                                                  ||"
    #echo "========================================================================"
    #
    #casa --nologger -c "flagdata(vis='${obsid}.ms',antenna='$ANT_FLAG',mode='manual',action='apply')"
    
    # Plotting the amplitude of the channels vs frequency.
    echo "========================================================================"
    echo "|| Plotting channels and baselines.                                   ||"
    echo "========================================================================"
    casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='frequency',yaxis='amp',correlation='xx,yy',ydatacolumn='data',coloraxis='spw',plotfile='amp-v-freq_${obsid}_data.png',showgui=False,overwrite=True)"
    casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='baseline',yaxis='amp',correlation='xx,yy',ydatacolumn='data',coloraxis='spw',plotfile='amp-v-base_${obsid}_data.png',showgui=False,overwrite=True)"

    # Deleting the gpu.fits files, these will just take up memory.
    rm -r *.fits

    # Return to the parent directory and restart the process for the next observation.
    cd ../../
done

exit 0
