
source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

input="mos1_target_ra.txt"
sky_ra=`cat $input`
input="mos1_target_dec.txt"
sky_dec=`cat $input`
esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset=mos1-fov-img.fits > mos1_esky2det.txt

input="mos1_target_ra.txt"
sky_ra=`cat $input`
input="mos1_target_dec_offset.txt"
sky_dec=`cat $input`
esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset=mos1-fov-img.fits > mos1_esky2det_offset1.txt

input="mos1_target_ra_offset.txt"
sky_ra=`cat $input`
input="mos1_target_dec.txt"
sky_dec=`cat $input`
esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset=mos1-fov-img.fits > mos1_esky2det_offset2.txt

input="mos2_target_ra.txt"
sky_ra=`cat $input`
input="mos2_target_dec.txt"
sky_dec=`cat $input`
esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset=mos2-fov-img.fits > mos2_esky2det.txt

input="mos2_target_ra.txt"
sky_ra=`cat $input`
input="mos2_target_dec_offset.txt"
sky_dec=`cat $input`
esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset=mos2-fov-img.fits > mos2_esky2det_offset1.txt

input="mos2_target_ra_offset.txt"
sky_ra=`cat $input`
input="mos2_target_dec.txt"
sky_dec=`cat $input`
esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset=mos2-fov-img.fits > mos2_esky2det_offset2.txt

