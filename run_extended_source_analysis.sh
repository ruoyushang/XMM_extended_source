
#source_name='extragalactic'
#ID_list=()
#ID_list+=('0827241101') # 51% SPF
#ID_list+=('0505460501') # 27% SPF
#ID_list+=('0804860301') # 25% SPF
#ID_list+=('0690900101') # 21% SPF
#ID_list+=('0803160301') # 10% SPF
#ID_list+=('0820560101') # 0% SPF
#ID_list+=('0827200401') # 4% SPF
#ID_list+=('0827200501') # 0% SPF
#ID_list+=('0827251001') # 0% SPF
#ID_list+=('0827251101') # 0% SPF


#source_name='Cas_A'
#ID_list=()
#ID_list+=('0137550301') # no MOS2
#ID_list+=('0400210101') # Cas A Northern lobe
#ID_list+=('0412180101') # Cas A
#ID_list+=('0782961401') # 34.7 arcmin away
#ID_list+=('0764640101') # 144.51 arcmin away, broken

#source_name='3HWC_J1928_p178'
#ID_list=()
#ID_list+=('0902120101')
#ID_list+=('0503740101')
#ID_list+=('0830192001')
#ID_list+=('0406730101')
#ID_list+=('0762980101')
#ID_list+=('0742710301')
#ID_list+=('0822330301')
#ID_list+=('0861880101') # 57% SPF
#ID_list+=('0841190101') # bright source
#ID_list+=('0851181701') # bright source
#ID_list+=('0724270101') # SNR with escaped CRs?
#ID_list+=('0724270201') # bright SNR
#ID_list+=('0742710401') # broken chip
#ID_list+=('0763240101') # strange ring patterns
#ID_list+=('0822330201') # strange ring patterns

#source_name='IC_443'
#ID_list=('ID0301960101')

source_name='skycoord_l0_b20'
ID_list=()
ID_list+=('0201770201')
ID_list+=('0651281401')
ID_list+=('0552070101')
ID_list+=('0763050301')
ID_list+=('0672660401')
ID_list+=('0852180401')
ID_list+=('0871390801')
ID_list+=('0800730901')
ID_list+=('0311590401')
ID_list+=('0784770801')
ID_list+=('0824420101')
ID_list+=('0502950101')
ID_list+=('0401870201')
ID_list+=('0094810301')
ID_list+=('0690330101')


for i in ${ID_list[@]}
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

