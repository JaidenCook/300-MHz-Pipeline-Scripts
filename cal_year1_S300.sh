#!/bin/bash -l

# Calibrate data using GLEAM year 1 catalogue

#SBATCH --account=pawsey0272
#SBATCH --partition=workq
#SBATCH --time=3:00:00
#SBATCH --nodes=1
##SBATCH --mem=60gb
#SBATCH --output=/astro/mwasci/jhue_cook/data/cal_year1.o%A_%a
#SBATCH --error=/astro/mwasci/jhue_cook/data/cal_year1.e%A_%a
#SBATCH --export=NONE
#SBATCH --array=1,21

start_time=`date +%s`

# Set aprun
#aprun="aprun -n 1 -d 20 -q "
#aprunsingle="aprun -n 1 -d 1 -q "
# now use and set srun
srun=" "
#srun -T 20
srunsinge=srun

# Set default values for optional parameters
flag=yes

# Read the options
TEMP=`getopt -o a:b:c:d:e:f --long input_data:,input_model:,output_dir:,obsid_list:,chan:,noflag -- "$@"`
eval set -- "$TEMP"

# Extract options and their arguments into variables
while true ; do
      case "$1" in
        -a|--input_data) # input directory containing raw measurement set (required argument)
            case "$2" in
                "") shift 2 ;;
                *) input_data=$2 ; shift 2 ;;
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
        -f|--noflag) flag=no ; shift ;; # do not flag tiles in obsid_list (no argument, acts as flag)
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# Check required arguments have been specified
if [ -z "$input_data" ]; then
  echo "Error: input_data not specified. Aborting."
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

# Set obsid
obsid=`sed "${SLURM_ARRAY_TASK_ID}q;d" $obsid_list | awk '{print $1}'`

# Set other input parameters
imsize=5000 # a quick image is made after calibration as a sanity check; imsize is the size of this image
#robust="-1.0"
robust="0.7"
ncpus=20

# Create output directory
if [ ! -e $output_dir ]; then
  mkdir $output_dir
fi

# Create snapshot directory
# remove 4 hashes below
if [ -e $output_dir/$obsid ]; then
  rm -rf $output_dir/$obsid
fi
mkdir $output_dir/$obsid
#cd $input_data/$obsid/
cd $output_dir/$obsid

# Write input parameters to file for record
cat >> input_parameters_cal_year1.txt <<EOPAR
input_data = $input_data
input_model = $input_model
output_dir = $output_dir
obsid_list = $obsid_list
chan = $chan
flag = $flag
imsize = $imsize
robust = $robust
ncpus = $ncpus
EOPAR

# Set pixel size
scale=`echo "0.55 / $chan" | bc -l`
scale=${scale:0:8}
#scale=0.04

# Copy measurement set, metafits file and sky model to output directory
if [ -e $input_data/$obsid/$obsid.ms ] && [ -e $input_data/$obsid/$obsid.metafits ] && [ -e $input_model/$obsid/skymodelformat.txt ]; then
# kluge to refresh stale file handles
    cd $input_data/$obsid/
    cd $output_dir/$obsid
    cp -r $input_data/$obsid/$obsid.ms $input_data/$obsid/$obsid.metafits $input_model/$obsid/skymodelformat.txt .
else
    echo "Error: input files are missing. Aborting."
    exit 1
fi

# -------------------------------------------------------------

# Flag tiles if required
if [ $flag == "yes" ]; then
  tile_list=`sed "${SLURM_ARRAY_TASK_ID}q;d" $obsid_list | awk '{print $2}'`
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

# determin max uv
maxuvm=887250/$chan

# Find calibration solutions
# edit out for now
echo "Running calibrate on  $obsid"
$srun calibrate -m skymodelformat.txt -minuv 60 -maxuv $maxuvm $obsid.ms ${obsid}_solutions.bin

# Apply the solutions
# edit out for now
echo "Applying solutions to  $obsid"
$srun applysolutions $obsid.ms ${obsid}_solutions.bin

# Really fast clean
echo "Doing a very fast clean (skipping"
# added Gausian tapering: 75"  (per J Cook)
#$srun wsclean -name ${obsid}_quick -size $imsize $imsize -niter 4000 -threshold 0.01 -pol xx,yy,xy,yx -weight briggs $robust -taper-gaussian 75 -scale $scale -stop-negative -small-inversion -join-polarizations -j $ncpus $obsid.ms
$srun wsclean -name ${obsid}_quick -size $imsize $imsize -niter 4000 -threshold 0.01 -pol xx,yy,xy,yx -weight briggs $robust -scale $scale -stop-negative -small-inversion -join-polarizations -j $ncpus $obsid.ms

# Plot phase and amplitude calibration solutions
echo "Doing plot of phase and amplityde"
$srunsingle aocal_plot.py --refant=127 ${obsid}_solutions.bin

# Re-plot amplitude calibration solutions, this time setting maximum of y axis to 100 and 10000
for amp in 100 10000; do
  mkdir t
  $srunsingle aocal_plot.py --refant=127 --outdir=./t --amp_max=$amp ${obsid}_solutions.bin
  mv ./t/${obsid}_solutions_amp.png ${obsid}_solutions_amp_max${amp}.png
  rm -rf t
done

# Make greyscale images
# no miriad installed so skip for now
#for pol in XX YY; do
## Extract rms noise
#  image=${obsid}_quick-${pol}-image.fits
#  rms=$(rms_measure.py -f $image)
#  pyhead.py -u IMAGERMS $rms $image
## Make greyscale maps using cgdisp in Miriad
#  fits in=$image out=map op=xyin
#  max=`echo "$rms*10.0" | bc -l`
#  min=`echo "$max*-1.0" | bc -l`
#  echo "max = $max"
#  echo "min = $min"
#  cgdisp device=/ps in=map labtyp=hms,dms options=blacklab,wedge range=$max,$min,lin,-1
## Convert map from ps to pdf format
#  ps2pdf pgplot.ps ${image%%.fits}.pdf
#  rm -rf map pgplot.ps
#done

# -------------------------------------------------------------

end_time=`date +%s`
duration=`echo "$end_time-$start_time" | bc -l`
echo "Total runtime = $duration sec"

# Move output and error files to output directory
mv $MYDATA/cal_year1.o${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} $MYDATA/cal_year1.e${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} .

exit 0
