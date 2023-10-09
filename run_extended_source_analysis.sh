
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

#source_name='skycoord_l0_b60'
#ID_list=()
#ID_list+=('0770580401')
#ID_list+=('0044740201')
#ID_list+=('0770580501')
#ID_list+=('0743900301')
#ID_list+=('0762610501')
#ID_list+=('0081340201')
#ID_list+=('0111281001')

source_name='PSR_J1928_p1746'
ID_list=()
ID_list+=('0764820101')
ID_list+=('0605370101')
ID_list+=('0742710301')
ID_list+=('0503740101')
ID_list+=('0554690101')
ID_list+=('0871190201')
ID_list+=('0504230201')
ID_list+=('0853980101')
ID_list+=('0841190101')
ID_list+=('0406730101')
ID_list+=('0830192001')
ID_list+=('0902120101')

#source_name='MGRO_J1908_p06'
#ID_list=()
#ID_list+=('0553640701')
#ID_list+=('0553640201')
#ID_list+=('0822330301')
#ID_list+=('0784040201')
#ID_list+=('0840910801')
#ID_list+=('0605700201')
#ID_list+=('0742820101')
#ID_list+=('0553640801')
#ID_list+=('0084100401')
#ID_list+=('0651680501')
#ID_list+=('0084100501')
#ID_list+=('0851181701')
#ID_list+=('0724270101')
#ID_list+=('0722250201')
#ID_list+=('0692210301')
#ID_list+=('0830450101')
#ID_list+=('0553640101')
#ID_list+=('0506430101')

#source_name='MGRO_J2019_p37'
#ID_list=()
#ID_list+=('0212481201')
#ID_list+=('0744640101')
#ID_list+=('0692810501')
#ID_list+=('0670480501')
#ID_list+=('0652770101')
#ID_list+=('0782961201')
#ID_list+=('0795712301')
#ID_list+=('0304000201')
#ID_list+=('0510011401')
#ID_list+=('0670480401')
#ID_list+=('0880030501')
#ID_list+=('0206240801')
#ID_list+=('0600030201')
#ID_list+=('0552350101')
#ID_list+=('0206240401')
#ID_list+=('0692810301')
#ID_list+=('0206240201')
#ID_list+=('0206240701')
#ID_list+=('0674050101')
#ID_list+=('0692810401')
#ID_list+=('0692810601')
#ID_list+=('0723310201')
#ID_list+=('0721570101')
#ID_list+=('0720600101')

mkdir /Users/rshang/xmm_analysis/output_plots/$source_name
for i in ${ID_list[@]}
do
    obs_ID='ID'$i
    rm -r /Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    mkdir /Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    outdir=/Users/rshang/xmm_analysis/output_plots/$source_name/$obs_ID
    (python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos1' 'fov') 2>&1 | tee $outdir/mos1_output.log 
    (python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos2' 'fov') 2>&1 | tee $outdir/mos2_output.log
    #python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos1' 'fov' 
    #python3 xmm_extended_source_analysis.py $source_name $obs_ID 'mos2' 'fov'
done

