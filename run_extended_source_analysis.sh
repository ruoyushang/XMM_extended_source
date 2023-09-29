#
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0827251101' # Percentage of flaring time: 0.0%
#
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0827251001' # Percentage of flaring time: 0.0% 
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0690900101' # Percentage of flaring time: 0.0% 
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0505460501' # Percentage of flaring time: 60%
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0820560101' # Percentage of flaring time: 48.8%
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0803160301' # Percentage of flaring time: 28.4%
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0827241101' # Percentage of flaring time: 15.2%
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0827200401' # Percentage of flaring time: 11.9%
#python3 xmm_extended_source_analysis.py 'extragalactic' 'ID0827200501' # Percentage of flaring time: 6%

#python3 xmm_extended_source_analysis.py 'RX_J1241.5_p3250' 'ID0056020901'
#python3 xmm_extended_source_analysis.py 'RX_J0256.5_p0006' 'ID0056020301'
#python3 xmm_extended_source_analysis.py 'RX_J1713.7_m3941' 'ID0093670301'
#python3 xmm_extended_source_analysis.py 'MGRO_J1908_p06' 'ID0553640201'
#python3 xmm_extended_source_analysis.py '3HWC_J1928_p178' 'ID0902120101'
#
#python3 xmm_extended_source_analysis.py 'Cas_A' 'ID0412180101'
#python3 xmm_extended_source_analysis.py 'Cas_A' 'ID0400210101' # Cas A Northern lobe
#python3 xmm_extended_source_analysis.py 'Cas_A' 'ID0782961401' # angular distance to Cas A: 34.7 arcmin
#python3 xmm_extended_source_analysis.py 'Cas_A' 'ID0764640101' # angular distance to Cas A: 144.51 arcmin, broken
#
#python3 xmm_extended_source_analysis.py 'PSR_J1856_p0245' 'ID0302970201'


#source_name='extragalactic'
#ID_array=()
#ID_array+=('ID0827241101') # 51% SPF
#ID_array+=('ID0505460501') # 27% SPF
#ID_array+=('ID0804860301') # 25% SPF
#ID_array+=('ID0690900101') # 21% SPF
#ID_array+=('ID0803160301') # 10% SPF
#ID_array+=('ID0820560101') # 0% SPF
#ID_array+=('ID0827200401') # 4% SPF
#ID_array+=('ID0827200501') # 0% SPF
#ID_array+=('ID0827251001') # 0% SPF
#ID_array+=('ID0827251101') # 0% SPF


#source_name='Cas_A'
#ID_array=()
#ID_array+=('ID0137550301') # no MOS2
#ID_array+=('ID0400210101') # Cas A Northern lobe
#ID_array+=('ID0412180101') # Cas A
#ID_array+=('ID0782961401') # 34.7 arcmin away
#ID_array+=('ID0764640101') # 144.51 arcmin away, broken

#source_name='3HWC_J1928_p178'
#ID_array=()
#ID_array+=('ID0902120101')
#ID_array+=('ID0503740101')
#ID_array+=('ID0830192001')
#ID_array+=('ID0406730101')
#ID_array+=('ID0762980101')
#ID_array+=('ID0742710301')
#ID_array+=('ID0822330301')
#ID_array+=('ID0861880101') # 57% SPF
#ID_array+=('ID0841190101') # bright source
#ID_array+=('ID0851181701') # bright source
#ID_array+=('ID0724270101') # SNR with escaped CRs?
#ID_array+=('ID0724270201') # bright SNR
#ID_array+=('ID0742710401') # broken chip
#ID_array+=('ID0763240101') # strange ring patterns
#ID_array+=('ID0822330201') # strange ring patterns

#source_name='IC_443'
#ID_array=('ID0301960101')

source_name='skycoord_l0_b25'
ID_array=()
ID_array+=('0083000101')
ID_array+=('0672660401')
ID_array+=('0844050101')
ID_array+=('0112880101')
ID_array+=('0556200301')
ID_array+=('0110930101')
ID_array+=('0691320601')
ID_array+=('0202940201')
ID_array+=('0822470101')


for i in ${ID_array[@]}
do
    obs_ID='ID'$i
    rm -r /Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    mkdir /Users/rshang/xmm_analysis/output_plots/$source_name
    mkdir /Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    outdir=/Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    (python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos1' 'fov') 2>&1 | tee $outdir/mos1_output.log 
    (python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos2' 'fov') 2>&1 | tee $outdir/mos2_output.log
    #python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos1' 'fov' 
    #python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos2' 'fov'
done

