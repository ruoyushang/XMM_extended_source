
source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

mos1_fits=$(ls P*M1*EVL*)
mos2_fits=$(ls P*M2*EVL*)
echo "All files:"
echo $mos1_fits
echo $mos2_fits

mos1_fits=$(ls P*M1S001MIEVL*)
mos2_fits=$(ls P*M2S002MIEVL*)
#mos1_fits=$(ls P*M1S006MIEVL*)
#mos2_fits=$(ls P*M2S007MIEVL*)
echo "Copy files:"
echo $mos1_fits
echo $mos2_fits
cp $mos1_fits mos1.fits
cp $mos2_fits mos2.fits

