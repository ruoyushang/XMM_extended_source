

source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

reg_x=73
reg_y=1515
reg_r=3600
bkg_r=7200

#em_std_cut='(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM'
#mos_src_region_cut='((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r'))'
##analysis_cut=$em_std_cut' && '$mos_src_region_cut
#analysis_cut=$mos_src_region_cut
#
## extract spectrum
#evselect table=mos1-fov-evt.fits energycolumn='PI' withfilteredset=yes filteredset=mos1-fov-src.fits keepfilteroutput=yes filtertype='expression' expression="$analysis_cut" withspectrumset=yes spectrumset=mos1-fov-pi.fits imageset=mos1-roi-img.fits spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999
#backscale spectrumset=mos1-fov-pi.fits
#rmfgen rmfset=mos1-src-rmf.fits spectrumset=mos1-fov-pi.fits
#arfgen arfset=mos1-src-arf.fits spectrumset=mos1-fov-pi.fits withrmfset=yes rmfset=mos1-src-rmf.fits
#
#evselect table=mos2-fov-evt.fits energycolumn='PI' withfilteredset=yes filteredset=mos2-fov-src.fits keepfilteroutput=yes filtertype='expression' expression="$analysis_cut" withspectrumset=yes spectrumset=mos2-fov-pi.fits imageset=mos2-roi-img.fits spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999
#backscale spectrumset=mos2-fov-pi.fits
#rmfgen rmfset=mos2-src-rmf.fits spectrumset=mos2-fov-pi.fits
#arfgen arfset=mos2-src-arf.fits spectrumset=mos2-fov-pi.fits withrmfset=yes rmfset=mos2-src-rmf.fits


em_std_cut='(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM'
mos_src_region_cut='((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$bkg_r')) && !((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r'))'
#analysis_cut=$em_std_cut' && '$mos_src_region_cut
analysis_cut=$mos_src_region_cut

# extract spectrum
evselect table=mos1-fov-evt.fits energycolumn='PI' withfilteredset=yes filteredset=mos1-bkg-src.fits keepfilteroutput=yes filtertype='expression' expression="$analysis_cut" withspectrumset=yes spectrumset=mos1-bkg-pi.fits imageset=mos1-bkg-img.fits spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999

evselect table=mos2-fov-evt.fits energycolumn='PI' withfilteredset=yes filteredset=mos2-bkg-src.fits keepfilteroutput=yes filtertype='expression' expression="$analysis_cut" withspectrumset=yes spectrumset=mos2-bkg-pi.fits imageset=mos2-bkg-img.fits spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999

