

source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/


em_std_cut='(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM'
mos_src_region_cut='((DETX,DETY) in ANNULUS(0, 0, 0, 7000))'
analysis_cut=$em_std_cut' && '$mos_src_region_cut

# extract spectrum
evselect table=mos1S001.fits energycolumn='PI' withfilteredset=yes filteredset=mos1S001-fov-src.fits keepfilteroutput=yes filtertype='expression' expression="$analysis_cut" withspectrumset=yes spectrumset=mos1S001-fov-pi.fits spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999
evselect table=mos2S002.fits energycolumn='PI' withfilteredset=yes filteredset=mos2S002-fov-src.fits keepfilteroutput=yes filtertype='expression' expression="$analysis_cut" withspectrumset=yes spectrumset=mos2S002-fov-pi.fits spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999


backscale spectrumset=mos1S001-fov-pi.fits
rmfgen rmfset=mos1_src_rmf.fits spectrumset=mos1S001-fov-pi.fits
arfgen arfset=mos1_src_arf.fits spectrumset=mos1S001-fov-pi.fits withrmfset=yes rmfset=mos1_src_rmf.fits

backscale spectrumset=mos2S002-fov-pi.fits
rmfgen rmfset=mos2_src_rmf.fits spectrumset=mos2S002-fov-pi.fits
arfgen arfset=mos2_src_arf.fits spectrumset=mos2S002-fov-pi.fits withrmfset=yes rmfset=mos2_src_rmf.fits
