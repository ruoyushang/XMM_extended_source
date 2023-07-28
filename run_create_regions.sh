
source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/


convregion mode=2 shape=CIRCLE inregion=pointing_coord.txt outregion=det_coord_mos1.txt imagefile=mos1-fov-img.fits
convregion mode=2 shape=CIRCLE inregion=pointing_coord.txt outregion=det_coord_mos2.txt imagefile=mos2-fov-img.fits

# The input file, reg_sky_coord.txt, contains:
# circle 207.12 26.59 4.0
# !circle 207.12 26.59 2.0
# The output data must then be edited into a region expression:
# &&(((DETX,DETY) IN circle(-2063,5786,4800))
# &&!((DETX,DETY) IN circle(-2063,5786,2400)))
