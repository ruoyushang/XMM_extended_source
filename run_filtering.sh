
# Soft proton (SP) filtering is accomplished with the espfilt task. i
# This task creates two light curves (one from the FOV data, and one from the corner data) in the 2.5-8.5 keV band, 
# and creates an X-ray count rate histogram from the FOV data. 
# For a typical observation the histogram will have a roughly Gaussian peak at some nominal count rate 
# (the count rate during time intervals unaffected, or at least minimally affected, by SP contamination) with a higher count-rate tail.
# espfilt will fit a reasonable Gaussian to the peak and determine thresholds at plus or minus 1.5σ. 
# The espfilt task then creates a GTI file for those time intervals with count rates within the thresholds 
# and uses the task evselect to filter the data to create “cleaned” photon event files.

source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

# Output files:
# mos1S001-allevc.fits
# mos2S002-allevc.fits
# pnS003-allevc.fits

rm mos1-*.fits
rm mos1-*.qdp
rm mos2-*.fits
rm mos2-*.qdp

espfilt eventfile=mos1.fits elow=$ELOW ehigh=$EHIGH withsmoothing=yes smooth=51 rangescale=6.0 allowsigma=1.5 method=histogram keepinterfiles=false
espfilt eventfile=mos2.fits elow=$ELOW ehigh=$EHIGH withsmoothing=yes smooth=51 rangescale=6.0 allowsigma=1.5 method=histogram keepinterfiles=false
mos1_gti=$(ls mos1*-gti.fits)
mos2_gti=$(ls mos2*-gti.fits)
mv $mos1_gti mos1-gti-tight.fits
mv $mos2_gti mos2-gti-tight.fits

espfilt eventfile=mos1.fits elow=$ELOW ehigh=$EHIGH withsmoothing=yes smooth=51 rangescale=6.0 allowsigma=2.5 method=histogram keepinterfiles=false
espfilt eventfile=mos2.fits elow=$ELOW ehigh=$EHIGH withsmoothing=yes smooth=51 rangescale=6.0 allowsigma=2.5 method=histogram keepinterfiles=false
mos1_gti=$(ls mos1*-gti.fits)
mos2_gti=$(ls mos2*-gti.fits)
mv $mos1_gti mos1-gti-medium.fits
mv $mos2_gti mos2-gti-medium.fits

espfilt eventfile=mos1.fits elow=$ELOW ehigh=$EHIGH withsmoothing=yes smooth=51 rangescale=6.0 allowsigma=3.5 method=histogram keepinterfiles=false
espfilt eventfile=mos2.fits elow=$ELOW ehigh=$EHIGH withsmoothing=yes smooth=51 rangescale=6.0 allowsigma=3.5 method=histogram keepinterfiles=false
mos1_gti=$(ls mos1*-gti.fits)
mos2_gti=$(ls mos2*-gti.fits)
mv $mos1_gti mos1-gti-loose.fits
mv $mos2_gti mos2-gti-loose.fits

#espfilt eventfile=pnS003.fits elow=$ELOW ehigh=$EHIGH withsmoothing=yes smooth=51 rangescale=15.0 allowsigma=3.0 method=histogram withoot=Y ootfile=pnS003-oot.fits keepinterfiles=false

# (FLAG & 0x766aa000)==0 is not equivalent to #XMMEA_EM but is the equivalent to #XMMEA EM with the addition of the out-of-FOV events
# (FLAG & 0x766ba000)==0 equivalent to #XMMEA_EM

# select FoV events
#em_fov_cut='(PATTERN<=12)&&(PI in ['$ELOW':'$EHIGH'])&&(FLAG & 0x766ba000)==0'
em_fov_cut='(PI in ['$ELOW':'$EHIGH'])&&(FLAG & 0x766ba000)==0'
evselect table=mos1.fits withfilteredset=yes filtertype=expression expression="$em_fov_cut" filteredset=mos1-fov-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2.fits withfilteredset=yes filtertype=expression expression="$em_fov_cut" filteredset=mos2-fov-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes


# select corner events
#em_cor_cut='(PATTERN<=12)&&(PI in ['$ELOW':'$EHIGH'])&&((FLAG & 0x766aa000)==0)&&((FLAG & 0x766ba000)!=0)'
em_cor_cut='(PI in ['$ELOW':'$EHIGH'])&&((FLAG & 0x766aa000)==0)&&((FLAG & 0x766ba000)!=0)'
evselect table=mos1.fits withfilteredset=yes filtertype=expression expression="$em_cor_cut" filteredset=mos1-cor-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes
evselect table=mos2.fits withfilteredset=yes filtertype=expression expression="$em_cor_cut" filteredset=mos2-cor-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes



# make the attitude file
atthkgen atthkset=attitude.fits timestep=1

# select cosmic-ray events
evqpb table=mos1.fits attfile=attitude.fits outset=mos1-fwc-evt.fits
evqpb table=mos2.fits attfile=attitude.fits outset=mos2-fwc-evt.fits

# select FoV events
em_fwc_fov_cut='(PI in ['$ELOW':'$EHIGH'])&&((FLAG & 0x10000)==0)'
evselect table=mos1-fwc-evt.fits withfilteredset=yes filtertype=expression expression="$em_fwc_fov_cut" filteredset=mos1-fwc-fov-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes
evselect table=mos2-fwc-evt.fits withfilteredset=yes filtertype=expression expression="$em_fwc_fov_cut" filteredset=mos2-fwc-fov-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes

# select corner events
em_fwc_cor_cut='(PI in ['$ELOW':'$EHIGH'])&&((FLAG & 0x10000)!=0)'
evselect table=mos1-fwc-evt.fits withfilteredset=yes filtertype=expression expression="$em_fwc_cor_cut" filteredset=mos1-fwc-cor-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes
evselect table=mos2-fwc-evt.fits withfilteredset=yes filtertype=expression expression="$em_fwc_cor_cut" filteredset=mos2-fwc-cor-evt.fits keepfilteroutput=yes updateexposure=yes filterexposure=yes



# to view the image
# imgdisplay withimagefile=true imagefile=mos1image_src_filt.fits
# dsplot x=TIME y=RATE table=mos1S001-fovlc.fits
