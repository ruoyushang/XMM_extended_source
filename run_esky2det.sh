
source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

echo $PWD
cd analysis/sky_det_ref/

for (( i=0; i <= 1; ++i ))
do
    for (( j=0; j <= 1; ++j ))
    do
        input="mos1_target_ra_${i}.txt"
        sky_ra=`cat $input`
        input="mos1_target_dec_${j}.txt"
        sky_dec=`cat $input`
        echo $sky_ra
        echo $sky_dec
        esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset='../mos1-fov-img.fits'
        esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset='../mos1-fov-img.fits' >"mos1_esky2det_ra_${i}_dec_${j}.txt"
    done
done

for (( i=0; i <= 1; ++i ))
do
    for (( j=0; j <= 1; ++j ))
    do
        input="mos2_target_ra_${i}.txt"
        sky_ra=`cat $input`
        input="mos2_target_dec_${j}.txt"
        sky_dec=`cat $input`
        echo $sky_ra
        echo $sky_dec
        esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset='../mos2-fov-img.fits'
        esky2det datastyle=user ra=$sky_ra dec=$sky_dec outunit=det calinfoset='../mos2-fov-img.fits' >"mos2_esky2det_ra_${i}_dec_${j}.txt"
    done
done

#esky2det datastyle=user ra=350.82 dec=58.82 outunit=det calinfoset='../mos1-fov-img.fits'
#esky2det datastyle=user ra=350.82 dec=58.82 outunit=det calinfoset='../mos1-fov-img.fits' >"mos1_esky2det_ra_97_dec_97.txt"
#esky2det datastyle=user ra=350.92 dec=58.82 outunit=det calinfoset='../mos1-fov-img.fits'
#esky2det datastyle=user ra=350.92 dec=58.82 outunit=det calinfoset='../mos1-fov-img.fits' >"mos1_esky2det_ra_97_dec_98.txt"
#esky2det datastyle=user ra=350.82 dec=58.92 outunit=det calinfoset='../mos1-fov-img.fits'
#esky2det datastyle=user ra=350.82 dec=58.92 outunit=det calinfoset='../mos1-fov-img.fits' >"mos1_esky2det_ra_98_dec_97.txt"
#
#esky2det datastyle=user ra=350.82 dec=58.82 outunit=det calinfoset='../mos2-fov-img.fits'
#esky2det datastyle=user ra=350.82 dec=58.82 outunit=det calinfoset='../mos2-fov-img.fits' >"mos2_esky2det_ra_97_dec_97.txt"
#esky2det datastyle=user ra=350.92 dec=58.82 outunit=det calinfoset='../mos2-fov-img.fits'
#esky2det datastyle=user ra=350.92 dec=58.82 outunit=det calinfoset='../mos2-fov-img.fits' >"mos2_esky2det_ra_97_dec_98.txt"
#esky2det datastyle=user ra=350.82 dec=58.92 outunit=det calinfoset='../mos2-fov-img.fits'
#esky2det datastyle=user ra=350.82 dec=58.92 outunit=det calinfoset='../mos2-fov-img.fits' >"mos2_esky2det_ra_98_dec_97.txt"


