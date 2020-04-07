#!/bin/bash

start_time=`date +%s`
# For 40kHz.
coarse_chans="0,32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,512,544,576,608,640,672,704,736,767"

# Read the options
TEMP=`getopt -o a:b:c: --long obsid:,data_column:,flag_chans: -- "$@"`
eval set -- "$TEMP"

# Extract options and their arguments into variables
while true ; do
      case "$1" in
        -a|--obsid) # input directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) obsid=$2 ; shift 2 ;;
            esac ;;
        -b|--data_column) # input directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) data_column=$2 ; shift 2 ;;
            esac ;;
        -c|--flag_chans) # input directory (required argument)
            case "$2" in
                "") shift 2 ;;
                *) flag_chans=$2 ; shift 2 ;;
            esac ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# Checking input obsid.
if [ -z "$obsid" ]; then
	echo "No obsid given exiting program."
	exit 0
fi

# Checking data_column input name.
if [ $data_column == "CORRECTED" ] || [ $data_column == "corrected" ]; then
	column_name="CORRECTED_DATA"
	echo "Flagging CORRECTED_DATA"
elif [ $data_column == "DATA" ] || [ $data_column == "data" ]; then
	column_name="DATA"
	echo "Flagging DATA"
else
	echo "Not a valid data column, should be CORRECTED or DATA"
	exit 0
fi

# Flagging specific coarse channels.
if [ -z "$flag_chans" ]; then
	echo "No coarse channels to flag."
else
	echo "======================================================================================="
	echo "Flagging coarse channel(s) [$flag_chans]."
	echo "======================================================================================="

	Nchans=$(echo $flag_chans | awk -F, '{print NF; exit}')
	echo "Number of channels to flag $Nchans"
	Nseq=$(seq $Nchans)

	for i in $Nseq;
	do
		chan_lower=$(echo $flag_chans | awk -v var="$i" -F',' '{ print $var }')
		chan_upper=$(("$chan_lower" + "1"))
		fine_chan_lower=$(echo $coarse_chans | awk -v var_lower=$chan_lower -F, '{ print $var_lower }')
		fine_chan_upper=$(echo $coarse_chans | awk -v var_upper=$chan_upper -F, '{ print $var_upper }')
		echo "#${i} ${fine_chan_lower} ${fine_chan_upper}"

		# Flagging specific channels
		casa --nologger -c "flagdata('${obsid}.ms',mode='clip',datacolumn='${data_column}',spw='0:${fine_chan_lower}~${fine_chan_upper}',clipminmax=[0.0,0.0],clipzeros=True,action='apply')"
	done
fi

# Further flagging of the corrected data is required to eliminate RFI.                                                                  
echo "======================================================================================="
echo "Baseline and RFI flagging."
echo "======================================================================================="

# Performing further RFI flagging.                                                                                                      
echo "Applying rflag."
echo "======================================================================================="

casa --nologger -c "flagdata('${obsid}.ms',mode='rflag',datacolumn='${data_column}',timedevscale=4.0,action='apply')"

echo "======================================================================================="
echo "Applying tcrop."
echo "======================================================================================="

casa --nologger -c "flagdata('${obsid}.ms',mode='tfcrop',datacolumn='${data_column}',freqcutoff=3.0,usewindowstats='sum',action='apply')"

echo "======================================================================================="
echo "Antennas to be flagged."
echo "======================================================================================="

ANT_FLAG=$(ms_flag_by_uvdist.py "${obsid}.ms" ${column_name})
echo "Antenna pairs to be flagged."
echo $ANT_FLAG

echo "======================================================================================="
echo "Applying steflag."
echo "======================================================================================="

# Flagging baselines.                                                                                                                   
casa --nologger -c "flagdata(vis='${obsid}.ms',antenna='$ANT_FLAG',mode='manual',action='apply')"

echo "======================================================================================="
echo "Plotting CORRECTED_DATA channels and baselines."
echo "======================================================================================="

casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='frequency',yaxis='amp',correlation='xx,yy',ydatacolumn='${data_column}',coloraxis='spw',plotfile='amp-v-freq_${obsid}_${data_column}-flgd.png',showgui=False,overwrite=True)"
casa --nogui --agg --nologger -c "plotms(vis='${obsid}.ms',xaxis='baseline',yaxis='amp',correlation='xx,yy',ydatacolumn='${data_column}',coloraxis='spw',plotfile='amp-v-base_${obsid}_${data_column}-flgd.png',showgui=False,overwrite=True)"

