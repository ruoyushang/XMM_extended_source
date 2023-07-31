
source set_env_after_cifbuild.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

cd analysis/

#src_x=73
#src_y=1515
#src_r=4000
#
#src_x=12066
#src_y=12905
#src_r=4000

src_x=0
src_y=0
src_r=4000

reg_x=$src_x
reg_y=$src_y
reg_r_inner=0
reg_r_outer=$src_r
region_cut='((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_outer')) && !((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_inner'))'

out_mos1_evt_file='mos1-reg1-evt.fits'
out_mos2_evt_file='mos2-reg1-evt.fits'
evselect table=mos1-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

out_mos1_evt_file='mos1-fwc-reg1-evt.fits'
out_mos2_evt_file='mos2-fwc-reg1-evt.fits'
evselect table=mos1-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

reg_r_inner=$src_r
reg_r_outer=$src_r+2000
region_cut='((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_outer')) && !((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_inner'))'

out_mos1_evt_file='mos1-reg2-evt.fits'
out_mos2_evt_file='mos2-reg2-evt.fits'
evselect table=mos1-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

out_mos1_evt_file='mos1-fwc-reg2-evt.fits'
out_mos2_evt_file='mos2-fwc-reg2-evt.fits'
evselect table=mos1-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

reg_r_inner=$src_r+2000
reg_r_outer=$src_r+4000
region_cut='((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_outer')) && !((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_inner'))'

out_mos1_evt_file='mos1-reg3-evt.fits'
out_mos2_evt_file='mos2-reg3-evt.fits'
evselect table=mos1-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

out_mos1_evt_file='mos1-fwc-reg3-evt.fits'
out_mos2_evt_file='mos2-fwc-reg3-evt.fits'
evselect table=mos1-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

reg_r_inner=$src_r+4000
reg_r_outer=$src_r+12000
region_cut='((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_outer')) && !((DETX,DETY) in CIRCLE('$reg_x','$reg_y','$reg_r_inner'))'

out_mos1_evt_file='mos1-reg4-evt.fits'
out_mos2_evt_file='mos2-reg4-evt.fits'
evselect table=mos1-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

out_mos1_evt_file='mos1-fwc-reg4-evt.fits'
out_mos2_evt_file='mos2-fwc-reg4-evt.fits'
evselect table=mos1-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos1_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes 
evselect table=mos2-fwc-fov-evt.fits withfilteredset=yes filtertype=expression expression="$region_cut" filteredset=$out_mos2_evt_file keepfilteroutput=yes updateexposure=yes filterexposure=yes

