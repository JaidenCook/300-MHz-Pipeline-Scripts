#!/bin/bash -l

# Build sky model for snapshot using GLEAM year 1 catalogue

#SBATCH --account=pawsey0272
#SBATCH --partition=workq
#SBATCH --time=7:00:00
#SBATCH --nodes=1
##SBATCH --mem=62gb
#SBATCH --output=/astro/mwasci/jhue_cook/data/build_beam.o%A_%a
#SBATCH --error=/astro/mwasci/jhue_cook/data/build_beam.e%A_%a
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
pbeam_size=5000
avg_yn="n"

# Read the options
TEMP=`getopt -o a:b:c:d:e:f:g: --long input_dir:,output_dir:,obsid_list:,chan:,pbeam_model:,pbeam_size:,avg_yn: -- "$@"`
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
        -e|--pbeam_model) # primary beam model (optional argument); set to 2014 or 2016
            case "$2" in
                "") shift 2 ;;
                *) pbeam_model=$2 ; shift 2 ;;
            esac ;;
        -f|--pbeam_size) # Image size for the primary beam, input is number of pixels.
            case "$2" in
                "") shift 2 ;;
                *) pbeam_size=$2 ; shift 2 ;;
            esac ;;
	-g|--avg_yn) # Image size for the primary beam, input is number of pixels.
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
elif [ -z "$output_dir" ]; then
  echo "Error: output_dir not specified. Aborting."
  exit 1
elif [ -z "$obsid_list" ]; then
  echo "Error: obsid_list not specified. Aborting."
  exit 1
elif [ -z "$chan" ]; then
  echo "Error: chan not specified. Aborting."
  exit 1
elif [ -z "$pbeam_size" ]; then
  echo "Error: pbeam_size not specified. Aborting."
  exit 1
fi

# Check chan parameter
if [ $chan != "69" ] && [ $chan != "93" ] && [ $chan != "121" ] && [ $chan != "145" ] && [ $chan != "169" ] && [ $chan != "236" ]; then
  echo "Error: chan should be set to 69, 93, 121, 145, 169 or 236. Aborting."
  exit 1
fi

# Check pbeam_model parameter
if [ $pbeam_model != "2014" ] && [ $pbeam_model != "2016" ]; then
  echo "Error: pbeam_model should be set to 2014 or 2016. Aborting."
  exit 1
fi

# Set obsid
obsid=`sed "${SLURM_ARRAY_TASK_ID}q;d" $obsid_list | awk '{print $1}'`

# Set other input parameters
absmem=64
ncpus=20

# Create output directory
if [ ! -e $output_dir ]; then
  mkdir $output_dir
fi

# Create snapshot directory
if [ -e $output_dir/$obsid ]; then
  echo "Error: Output directory $output_dir/$obsid already exists. Aborting."
  exit 1
else
  mkdir $output_dir/$obsid
  cd $output_dir/$obsid
fi

# Write input parameters to file for record
cat >> input_parameters_build_beam_S300.txt <<EOPAR
input_dir = $input_dir
output_dir = $output_dir
obsid_list = $obsid_list
chan = $chan
pbeam_model = $pbeam_model
pbeam_size = $pbeam_size
absmem = $absmem
ncpus = $ncpus
avg_yn = $avg_yn
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

# Set pixel size
scale=`echo "8.25 / $chan" | bc -l`
scale=${scale:0:8}
#scale=0.04

#I am getting errors here, CHANGE IMMEDIATELY.
# Copy measurement set and metafits file to output directory
if [ -e $input_dir/$obsid/$obsid.ms ] && [ -e $input_dir/$obsid/$obsid.metafits ]; then
   ln -s $input_dir/$obsid/$obsid.ms $input_dir/$obsid/$obsid.metafits .
else
  echo "Error: input files are missing. Aborting."
  exit 1
fi

# -------------------------------------------------------------
# Make PSF image at different time intervals.

# Make PSF image (need this as a template for creating primary beam)
if [ $avg_yn == "y" ]; then
	echo "Making time averaged beam"
	echo "Making initial time sampled PB."
	echo "$srun wsclean -name $obsid -size $pbeam_size $pbeam_size -weight briggs -1.0 -scale $scale -absmem $absmem -smallinversion -j $ncpus -make-psf-only -interval 1 3 $obsid.ms"

	$srun wsclean -name $obsid -size $pbeam_size $pbeam_size -weight briggs 0.7 -scale $scale -abs-mem $absmem -small-inversion -j $ncpus -make-psf-only -interval 1 3 $obsid.ms

	echo "Making Beam"

	$srun make_beam.py -v -f ${obsid}-psf.fits -m $obsid.metafits --model=$pbeam_model

	# Make Stokes I primary beam model by taking mean of Stokes XX & YY primary beam models
	$srun calc_stokesi_beam.py<<-EOF
	${obsid}-psf_beamXX.fits
	${obsid}-psf_beamYY.fits
	${obsid}-psf_beamI.fits
	EOF

	rm ${obsid}-psf_beamXX.fits ${obsid}-psf_beamYY.fits

	mv ${obsid}-psf_beamI.fits ${obsid}-psf_beamI_ti.fits

	echo "#-----------------------------------------------------------"
	echo "Making final time sampled PB."
	echo "$srun wsclean -name $obsid -size $pbeam_size $pbeam_size -weight briggs -1.0 -scale $scale -absmem $absmem -smallinversion -j $ncpus -make-psf-only -interval 211 213 $obsid.ms"
	
	$srun wsclean -name $obsid -size $pbeam_size $pbeam_size -weight briggs 0.7 -scale $scale -abs-mem $absmem -small-inversion -j $ncpus -make-psf-only -interval 211 213 $obsid.ms
	
	echo "Making Beam"
	
	$srun make_beam.py -v -f ${obsid}-psf.fits -m $obsid.metafits --model=$pbeam_model
	
	
	# Make Stokes I primary beam model by taking mean of Stokes XX & YY primary beam models
	$srun calc_stokesi_beam.py<<-EOF
	${obsid}-psf_beamXX.fits
	${obsid}-psf_beamYY.fits
	${obsid}-psf_beamI.fits
	EOF

	rm ${obsid}-psf_beamXX.fits ${obsid}-psf_beamYY.fits

	mv ${obsid}-psf_beamI.fits ${obsid}-psf_beamI_tf.fits

elif [ $avg_yn == "n" ]; then
	echo "Not building average beam."
fi

echo "#-----------------------------------------------------------"
echo "Building the sub-channel beams"
# Make PSF image (need this as a template for creating primary beam)
echo "$srun wsclean -name $obsid -size $pbeam_size $pbeam_size -weight briggs -1.0 -scale $scale -absmem $absmem -smallinversion -j $ncpus -make-psf-only -channelsout 4 $obsid.ms"

#$srun wsclean -name $obsid -size $pbeam_size $pbeam_size -weight briggs 0.7 -scale $scale -abs-mem $absmem -small-inversion -j $ncpus -make-psf-only -channels-out 8 $obsid.ms

$srun wsclean -name $obsid -size $pbeam_size $pbeam_size -weight briggs 0.7 -scale $scale -abs-mem $absmem -small-inversion -j $ncpus -make-psf-only -channels-out 4 $obsid.ms

# Make Stokes XX & YY primary beam models for MFS, 0000 and 0003 subbands
# Ammended this to return the primary beam models for all subbands
#for subband in MFS 0000 0001 0002 0003 0004 0005 0006 0007; do

#for subband in MFS; do
for subband in MFS 0000 0001 0002 0003; do

$srun make_beam.py -v -f ${obsid}-${subband}-psf.fits -m $obsid.metafits --model=$pbeam_model

rm ${obsid}-${subband}-psf.fits

# Make Stokes I primary beam model by taking mean of Stokes XX & YY primary beam models
$srun calc_stokesi_beam.py<<EOF
${obsid}-${subband}-psf_beamXX.fits
${obsid}-${subband}-psf_beamYY.fits
${obsid}-${subband}-psf_beamI.fits
EOF

#rm ${obsid}-${subband}-psf_beamXX.fits ${obsid}-${subband}-psf_beamYY.fits
done


if [ $avg_yn == "y" ]; then
	echo "Creating time averaged PB"
	$srun Make_Beam_avg.py --pbeam_ti ${obsid}-psf_beamI_ti.fits --pbeam ${obsid}-MFS-psf_beamI.fits --pbeam_tf ${obsid}-psf_beamI_tf.fits
elif [ $avg_yn == "n" ]; then
	echo "#-----------------------------------------------------------" 
fi

end_time=`date +%s`
duration=`echo "$end_time-$start_time" | bc -l`
echo "Total runtime = $duration sec"

# Move output and error files to output directory
mv $MYDATA/build_beam.o${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} $MYDATA/build_beam.e${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} .

exit 0

