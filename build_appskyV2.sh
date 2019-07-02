#!/bin/bash -l

# Build sky model for snapshot using GLEAM year 1 catalogue

#SBATCH --account=pawsey0272
#SBATCH --partition=workq
#SBATCH --time=3:00:00
#SBATCH --nodes=1
##SBATCH --mem=62gb
#SBATCH --output=/astro/mwasci/jhue_cook/data/build_appsky.o%A_%a
#SBATCH --error=/astro/mwasci/jhue_cook/data/build_appsky.e%A_%a
#SBATCH --export=NONE
#SBATCH --array=1-84

start_time=`date +%s`

# Set aprun
#aprun="aprun -n 1 -d 20 -q "
#aprunsingle="aprun -n 1 -d 1 -q "
# Set srun
srun=" "
srunsingle="srun -T 1"

#Here I am placing the stuff needed to run jobs on galaxy:
# --account=mwasci
# --account=pawsey0272

# Set default values for optional parameters
pbeam_model=2016
# reduce to find many sources 
threshold=1
resolution=1.2
avg_yn="n"

# Read the options
TEMP=`getopt -o a:b:c:d:e:f: --long input_dir:,output_dir:,obsid_list:,threshold:,chan:,avg_yn: -- "$@"`
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
	-d|--threshold) # obsID list (optional argument)
            case "$2" in
                "") shift 2 ;;
                *) threshold=$2 ; shift 2 ;;
            esac ;;
        -e|--chan) # channel (required argument); set to 69, 93, 121, 145 or 169
            case "$2" in
                "") shift 2 ;;
                *) chan=$2 ; shift 2 ;;
            esac ;;
	-f|--avg_yn) # Image size for the primary beam, input is number of pixels.
            case "$2" in
                "") shift 2 ;;
                *) avg_yn=$2 ; shift 2 ;;
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
if [ $chan != "69" ] && [ $chan != "93" ] && [ $chan != "121" ] && [ $chan != "145" ] && [ $chan != "169" ] && [ $chan != "236" ]; then
  echo "Error: chan should be set to 69, 93, 121, 145, 169 or 236. Aborting."
  exit 1
fi

# Set obsid
obsid=`sed "${SLURM_ARRAY_TASK_ID}q;d" $obsid_list | awk '{print $1}'`

# Set other input parameters
absmem=64
ncpus=20

#catalogue=/group/mwasci/tfranzen/GLEAM_IDR6/GLEAMIDR6_published.fits # to be used to calibrate the data

catalogue=/home/jhue_cook/S300_cat/Total_S300MHz_Cat.fits #Change 6: Catalogue changed to my catalogue.

# Create snapshot directory
if [ -e $output_dir/$obsid ]; then
  cd $output_dir/$obsid
#  echo "Error: Output directory $output_dir/$obsid already exists. Aborting."
#  exit 1
else
  cd $output_dir/$obsid
fi

# Write input parameters to file for record
cat >> input_parameters_build_appsky.txt <<EOPAR
input_dir = $input_dir
output_dir = $output_dir
obsid_list = $obsid_list
threshold = $threshold
resolution = $resolution
pbeam_size = $pbeam_size
absmem = $absmem
ncpus = $ncpus
catalogue = $catalogue
EOPAR

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
fi
echo "Central frequency for channel $chan = $freq MHz"

# Get RA and Dec of pointing centre from .metafits file, in deg
ra_pnt=$(fitshdr $obsid.metafits | grep 'RA of pointing center' | awk '{print $3}')
dec_pnt=$(fitshdr $obsid.metafits | grep 'Dec of pointing center' | awk '{print $3}')

echo $obsid

if [ $avg_yn == "y" ]; then
	echo "Using time averaged beam"
	# Run attenuate by beam.
	$srun attenuate_by_beam_S300V6.py --pbeam ${obsid}-MFS-psf_beamI_avg.fits --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold	
elif [ $avg_yn == "n" ]; then
	# Run attenuate by beam.
        echo "Not using time averaged beam"
	$srun attenuate_by_beam_S300V6.py --pbeam ${obsid}-MFS-psf_beamI.fits --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold
fi

# Convert the edited catalogue to a format which is readable by calibrate
# Sources with int_flux_wide/peak_flux_wide < resolution will be considered to be unresolved
$srun vo2newmodel_S300V2.py --catalogue model_morecolumns.vot --output skymodelformat.txt --freq $freq --fluxcol S_centralfreq_uncorrected --alphacol alpha_uncorrected --qcol q_uncorrected --point --resolution=$resolution

# Remove .ms and .metafits files and PSF template images
#rm -rf $obsid.ms $obsid.metafits *.fits

# -------------------------------------------------------------

end_time=`date +%s`
duration=`echo "$end_time-$start_time" | bc -l`
echo "Total runtime = $duration sec"

# Move output and error files to output directory
mv $MYDATA/build_appsky.o${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} $MYDATA/build_appsky.e${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} .

exit 0

