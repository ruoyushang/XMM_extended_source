
source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

# Output files:
# mos1S001.fits
# mos2S002.fits
# pnS003.fits
# pnS003-oot.fits


# The emchain task will produce event files for all CCDS for MOS1 and MOS2 for imaging exposures, i
# as well as the outer CCDs, 2 through 7, for timing observations.

emchain

# However, epchain, by default, will create event files for only the first imaging exposure
# and must be run explicitly calling out the exposure number if there are multiple exposures (e.g., epchain exposure=2)
# In the ODF directory, look for files of the form * PN*IME*.FIT. There will be one file for each of the 12 CCDs for each segment. 
# In the file name, the segment number (like S001) follows directly after the ‘PN’. 
# An even easier way to determine how many pn segments there are is to type “epchain exposure=99”; you will get a list of pn segments
#
#epchain exposure=99

#epchain
#pn_fits=$(ls P*PN*PIEVL*)
#echo $pn_fits
#cp $pn_fits pnS003.fits
#
#epchain withoutoftime=true
#pn_oot_fits=$(ls P*PN*OOEVL*)
#echo $pn_oot_fits
#cp $pn_oot_fits pnS003-oot.fits

#rm *.FIT
#rm *.dat
#rm *.img
#rm errors

