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
#python3 xmm_extended_source_analysis.py 'Cas_A' 'ID0764640101' # angular distance to Cas A: 144.51 arcmin
#
#python3 xmm_extended_source_analysis.py 'PSR_J1856_p0245' 'ID0302970201'


#source_name='extragalactic'
#ID_array=('ID0827251001')  # SP-free
#ID_array=('ID0827251101')
#ID_array=('ID0505460501' 'ID0690900101' 'ID0803160301' 'ID0820560101' 'ID0827200401' 'ID0827200501' 'ID0827241101' 'ID0827251001' 'ID0827251101')
#
#source_name='Cas_A'
#ID_array=('ID0400210101')
#ID_array=('ID0412180101' 'ID0400210101' 'ID0782961401' 'ID0764640101') # ON, Cas A Northern lobe, angular distance to Cas A: 34.7 arcmin, 144.51 arcmin

#source_name='3HWC_J1928_p178'
#ID_array=('ID0841190101')
#ID_array=('ID0902120101' 'ID0503740101' 'ID0830192001' 'ID0841190101' 'ID0861880101')

source_name='IC_443'
ID_array=('ID0301960101')


for i in ${ID_array[@]}
do
    obs_ID=$i
    rm -r /Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    mkdir /Users/rshang/xmm_analysis/output_plots/$source_name
    mkdir /Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    outdir=/Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    (python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos1' 'fov') 2>&1 | tee $outdir/mos1_output.log 
    (python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos2' 'fov') 2>&1 | tee $outdir/mos2_output.log
    #python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos1' 'fov' 
    #python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos2' 'fov'
done

