#!/bin/bash -l

echo "You rang!"

obsid=1126709848
threshold=0.0127
#catalogue=Total_S300MHz_Cat.fits
catalogue=Total_S300MHz_Cat.fits
#freq=300
#freq=204
#Might need to change resolution to be that of the 300MHz survey.
resolution=1.2

#echo $obsid

# Get RA and Dec of pointing centre from .metafits file, in deg
ra_pnt=$(fitshdr $obsid.metafits | grep 'RA of pointing center' | awk '{print $3}')
dec_pnt=$(fitshdr $obsid.metafits | grep 'Dec of pointing center' | awk '{print $3}')
freq=$(fitshdr $obsid.metafits | grep 'FREQCENT' | awk '{print $2}')
delays=$(fitshdr $obsid.metafits | grep 'DELAYS' | awk '{print $3}')

#echo "freq = $freq"
#echo "ra_pnt = $ra_pnt"
#echo "dec_pnt = $dec_pnt"
#echo "delays = $del"

# Run attenuate_by_beam.py. This will take the catalogue and select all sources which lie within a distance of 
# python S_central_uncor_test_V6_2gleam.py --pbeam ${obsid}-MFS-psf_beamI.fits --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold

# Run attenuate_by_beam.py. This will take the catalogue and select all sources which lie within a distance of 
#python Model_format_0.py --obsid $obsid --freq $freq --delays $delays --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold
#python Model_format_1.py --obsid $obsid --freq $freq --delays $delays --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold
#python Model_format_2.py --obsid $obsid --freq $freq --delays $delays --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold
python Model_format_3.py --obsid $obsid --freq $freq --delays $delays --catalogue $catalogue --ra $ra_pnt --dec $dec_pnt --threshold $threshold


# Convert the edited catalogue to a format which is readable by calibrate
# Sources with int_flux_wide/peak_flux_wide < resolution will be considered to be unresolved
#python vo2newmodel_S300V2.py --catalogue model_morecolumns_temp.vot --output skymodelformat.txt --freq $freq --fluxcol S_centralfreq --alphacol alpha --qcol q_curve --point --resolution=$resolution
#python vo2newmodel_S300V2.py --catalogue model_morecolumns_temp.vot --output skymodelformat.txt --freq $freq --fluxcol S_centralfreq_uncorrected --alphacol alpha_uncorrected --qcol q_uncorrected --point --resolution=$resolution
#python vo2newmodel_2.py --catalogue model_morecolumns_temp.vot --output skymodelformat.txt --freq $freq --fluxcol S_centralfreq_uncorrected --coeff apparent_poly_coeff --point --resolution=$resolution
python vo2newmodel_2.py --catalogue model_morecolumns_temp.vot --output skymodelformat.txt --freq $freq --fluxcol S_centralfreq_uncorrected --coeff apparent_poly_coeff --point --resolution=$resolution

# -------------------------------------------------------------
# Move output and error files to output directory
##mv $MYDATA/build_model.o${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} $MYDATA/build_model.e${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} .

exit 0

