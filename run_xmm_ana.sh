
#cd /Users/rshang/xmm_analysis/extragalactic/ID0827241101
#cd /Users/rshang/xmm_analysis/extragalactic/ID0505460501
#cd /Users/rshang/xmm_analysis/extragalactic/ID0690900101
#cd /Users/rshang/xmm_analysis/extragalactic/ID0803160301
#cd /Users/rshang/xmm_analysis/extragalactic/ID0820560101
#cd /Users/rshang/xmm_analysis/extragalactic/ID0827200401
#cd /Users/rshang/xmm_analysis/extragalactic/ID0827200501
#cd /Users/rshang/xmm_analysis/extragalactic/ID0827251001
#cd /Users/rshang/xmm_analysis/extragalactic/ID0827251101
#cd /Users/rshang/xmm_analysis/extragalactic/ID0804860301
#
#cd /Users/rshang/xmm_analysis/MGRO_J1908_p06/ID0553640201
#cd /Users/rshang/xmm_analysis/RX_J1241.5_p3250/ID0056020901
#cd /Users/rshang/xmm_analysis/RX_J0256.5_p0006/ID0056020301
#cd /Users/rshang/xmm_analysis/RX_J1713.7_m3941/ID0093670301
#
#cd /Users/rshang/xmm_analysis/Cas_A/ID0137550301
#cd /Users/rshang/xmm_analysis/Cas_A/ID0412180101
#cd /Users/rshang/xmm_analysis/Cas_A/ID0400210101
#cd /Users/rshang/xmm_analysis/Cas_A/ID0764640101 # broken
#cd /Users/rshang/xmm_analysis/Cas_A/ID0782961401
#
#cd /Users/rshang/xmm_analysis/IC_443/ID0301960101

#cd /Users/rshang/xmm_analysis/PSR_J1856_p0245/ID0302970201
#cd /Users/rshang/xmm_analysis/PSR_J1856_p0245/ID0505920101

#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0902120101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0503740101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0830192001
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0841190101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0861880101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0406730101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0762980101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0851181701
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0724270101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0724270201
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0742710301
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0742710401
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0763240101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0822330201
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0822330301

#source_name='MGRO_J1908_p06'
#obs_ID_list=()
#obs_ID_list+=('0553640701')
#obs_ID_list+=('0553640201')
#obs_ID_list+=('0822330301')
#obs_ID_list+=('0784040201')
#obs_ID_list+=('0840910801')
#obs_ID_list+=('0605700201')
#obs_ID_list+=('0742820101')
#obs_ID_list+=('0553640801')
#obs_ID_list+=('0084100401')
#obs_ID_list+=('0651680501')
#obs_ID_list+=('0084100501')
#obs_ID_list+=('0851181701')
#obs_ID_list+=('0724270101')
#obs_ID_list+=('0722250201')
#obs_ID_list+=('0692210301')
#obs_ID_list+=('0830450101')
#obs_ID_list+=('0553640101')
#obs_ID_list+=('0506430101')

source_name='MGRO_J2019_p37'
obs_ID_list=()
#obs_ID_list+=('0212481201')
obs_ID_list+=('0744640101')
#obs_ID_list+=('0692810501')
#obs_ID_list+=('0670480501')
#obs_ID_list+=('0652770101')
#obs_ID_list+=('0782961201')
#obs_ID_list+=('0795712301')
#obs_ID_list+=('0304000201')
#obs_ID_list+=('0510011401')
#obs_ID_list+=('0670480401')
#obs_ID_list+=('0880030501')
#obs_ID_list+=('0206240801')
#obs_ID_list+=('0600030201')
#obs_ID_list+=('0552350101')
#obs_ID_list+=('0206240401')
#obs_ID_list+=('0692810301')
#obs_ID_list+=('0206240201')
#obs_ID_list+=('0206240701')
#obs_ID_list+=('0674050101')
#obs_ID_list+=('0692810401')
#obs_ID_list+=('0692810601')
#obs_ID_list+=('0723310201')
#obs_ID_list+=('0721570101')
#obs_ID_list+=('0720600101')

for i in ${obs_ID_list[@]}
do
    obs_ID=$i
    obs_folder='ID'$obs_ID
    echo $obs_folder

    cd /Users/rshang/xmm_analysis/$source_name/$obs_folder
    cp /Users/rshang/xmm_analysis/analysis_code/*.sh .
    cp /Users/rshang/xmm_analysis/analysis_code/*.py .

    #sh extract_data.sh
    #sh run_cifbuild.sh
    #sh clean.sh 
    #sh run_ana_chains.sh
    
    sh run_rename.sh $1 $2
    
    #sh run_filtering.sh
    #sh run_skyref.sh
    #sh run_esky2det.sh
    #sh run_response.sh

done











#sh run_regions.sh

#python3 run_pointing_radec.py
#sh run_create_regions.sh
#
#sh run_src_detect.sh

#sh run_create_regions.sh
#python3 convert_regions.py
#
#sh run_extract_spec_from_region.sh

#sh run_specgroup.sh
#sh run_img_prod.sh

