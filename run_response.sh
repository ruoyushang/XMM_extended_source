

source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

out_mos1_evt_file='mos1-reg1-evt.fits'
out_mos2_evt_file='mos2-reg1-evt.fits'
#mosspectra eventfile=$out_mos1_evt_file keepinterfiles=yes pattern=12 elow=200 ehigh=12000 ccds="T T T F F F T"
#mosspectra eventfile=$out_mos2_evt_file keepinterfiles=yes pattern=12 elow=200 ehigh=12000 ccds="T T T T T T T"
mosback inspecfile=mos1S001-fovt.pi elow=200 ehigh=12000 ccds="T T T F F F T" 
mosback inspecfile=mos2S002-fovt.pi elow=200 ehigh=12000 ccds="T T T T T T T" 

