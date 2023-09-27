

source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

evselect table=mos1-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,0,4000))' withspectrumset=yes spectrumset=mos1-fov-pi.fits specchannelmin=0 specchannelmax=11999 
evselect table=mos2-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,0,4000))' withspectrumset=yes spectrumset=mos2-fov-pi.fits specchannelmin=0 specchannelmax=11999

rmfgen rmfset=mos1-fov-rmf.fits spectrumset=mos1-fov-pi.fits
arfgen arfset=mos1-fov-arf.fits spectrumset=mos1-fov-pi.fits withrmfset=yes rmfset=mos1-fov-rmf.fits extendedsource=yes
rmfgen rmfset=mos2-fov-rmf.fits spectrumset=mos2-fov-pi.fits
arfgen arfset=mos2-fov-arf.fits spectrumset=mos2-fov-pi.fits withrmfset=yes rmfset=mos2-fov-rmf.fits extendedsource=yes
