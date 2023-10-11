

source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

evselect table=mos1-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,0,1000))' withspectrumset=yes spectrumset=mos1-fov-r0-pi.fits specchannelmin=0 specchannelmax=11999 
evselect table=mos2-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,0,1000))' withspectrumset=yes spectrumset=mos2-fov-r0-pi.fits specchannelmin=0 specchannelmax=11999
#backscale spectrumset=mos1-fov-r0-pi.fits badpixlocation=mos1-fov-evt.fits
rmfgen rmfset=mos1-fov-r0-rmf.fits spectrumset=mos1-fov-r0-pi.fits
arfgen arfset=mos1-fov-r0-arf.fits spectrumset=mos1-fov-r0-pi.fits withrmfset=yes rmfset=mos1-fov-r0-rmf.fits extendedsource=yes
#backscale spectrumset=mos2-fov-r0-pi.fits badpixlocation=mos2-fov-evt.fits
rmfgen rmfset=mos2-fov-r0-rmf.fits spectrumset=mos2-fov-r0-pi.fits
arfgen arfset=mos2-fov-r0-arf.fits spectrumset=mos2-fov-r0-pi.fits withrmfset=yes rmfset=mos2-fov-r0-rmf.fits extendedsource=yes

evselect table=mos1-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,5000,1000))' withspectrumset=yes spectrumset=mos1-fov-r1-pi.fits specchannelmin=0 specchannelmax=11999 
evselect table=mos2-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,5000,1000))' withspectrumset=yes spectrumset=mos2-fov-r1-pi.fits specchannelmin=0 specchannelmax=11999
#backscale spectrumset=mos1-fov-r1-pi.fits badpixlocation=mos1-fov-evt.fits
rmfgen rmfset=mos1-fov-r1-rmf.fits spectrumset=mos1-fov-r1-pi.fits
arfgen arfset=mos1-fov-r1-arf.fits spectrumset=mos1-fov-r1-pi.fits withrmfset=yes rmfset=mos1-fov-r1-rmf.fits extendedsource=yes
#backscale spectrumset=mos2-fov-r1-pi.fits badpixlocation=mos2-fov-evt.fits
rmfgen rmfset=mos2-fov-r1-rmf.fits spectrumset=mos2-fov-r1-pi.fits
arfgen arfset=mos2-fov-r1-arf.fits spectrumset=mos2-fov-r1-pi.fits withrmfset=yes rmfset=mos2-fov-r1-rmf.fits extendedsource=yes

evselect table=mos1-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,10000,1000))' withspectrumset=yes spectrumset=mos1-fov-r2-pi.fits specchannelmin=0 specchannelmax=11999 
evselect table=mos2-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,10000,1000))' withspectrumset=yes spectrumset=mos2-fov-r2-pi.fits specchannelmin=0 specchannelmax=11999
#backscale spectrumset=mos1-fov-r2-pi.fits badpixlocation=mos1-fov-evt.fits
rmfgen rmfset=mos1-fov-r2-rmf.fits spectrumset=mos1-fov-r2-pi.fits
arfgen arfset=mos1-fov-r2-arf.fits spectrumset=mos1-fov-r2-pi.fits withrmfset=yes rmfset=mos1-fov-r2-rmf.fits extendedsource=yes
#backscale spectrumset=mos2-fov-r2-pi.fits badpixlocation=mos2-fov-evt.fits
rmfgen rmfset=mos2-fov-r2-rmf.fits spectrumset=mos2-fov-r2-pi.fits
arfgen arfset=mos2-fov-r2-arf.fits spectrumset=mos2-fov-r2-pi.fits withrmfset=yes rmfset=mos2-fov-r2-rmf.fits extendedsource=yes

evselect table=mos1-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,15000,1000))' withspectrumset=yes spectrumset=mos1-fov-r3-pi.fits specchannelmin=0 specchannelmax=11999 
evselect table=mos2-fov-evt.fits expression='((DETX,DETY) in CIRCLE(0,15000,1000))' withspectrumset=yes spectrumset=mos2-fov-r3-pi.fits specchannelmin=0 specchannelmax=11999
#backscale spectrumset=mos1-fov-r3-pi.fits badpixlocation=mos1-fov-evt.fits
rmfgen rmfset=mos1-fov-r3-rmf.fits spectrumset=mos1-fov-r3-pi.fits
arfgen arfset=mos1-fov-r3-arf.fits spectrumset=mos1-fov-r3-pi.fits withrmfset=yes rmfset=mos1-fov-r3-rmf.fits extendedsource=yes
#backscale spectrumset=mos2-fov-r3-pi.fits badpixlocation=mos2-fov-evt.fits
rmfgen rmfset=mos2-fov-r3-rmf.fits spectrumset=mos2-fov-r3-pi.fits
arfgen arfset=mos2-fov-r3-arf.fits spectrumset=mos2-fov-r3-pi.fits withrmfset=yes rmfset=mos2-fov-r3-rmf.fits extendedsource=yes

