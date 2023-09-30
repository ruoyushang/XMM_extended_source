
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

select_mos1_fits=$1
select_mos2_fits=$2

if [ -z "$select_mos1_fits" ]
then
      echo "\$select_mos1_fits is empty"
      select_mos1_fits=$(ls P*M1S001MIEVL*)
      select_mos2_fits=$(ls P*M2S002MIEVL*)
      echo "Copy default files:"
      echo $select_mos1_fits
      echo $select_mos2_fits
else
      echo "\$select_mos1_fits is NOT empty"
      echo "Copy input files:"
      echo $select_mos1_fits
      echo $select_mos2_fits
fi

#mos1_fits=$(ls P*M1S001MIEVL*)
#mos2_fits=$(ls P*M2S002MIEVL*)
#mos1_fits=$(ls P*M1S001MIEVL*)
#mos2_fits=$(ls P*M2S004MIEVL*)
#echo "Copy files:"
#echo $mos1_fits
#echo $mos2_fits
cp $select_mos1_fits mos1.fits
cp $select_mos2_fits mos2.fits

