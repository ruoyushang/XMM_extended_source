
cd /Users/rshang/xmm_analysis/extragalactic/ID0827251001
#
#cd /Users/rshang/xmm_analysis/RX_J1241.5_p3250/ID0056020901
#cd /Users/rshang/xmm_analysis/RX_J0256.5_p0006/ID0056020301
#cd /Users/rshang/xmm_analysis/RX_J1713.7_m3941/ID0093670301
#cd /Users/rshang/xmm_analysis/Cas_A/ID0412180101
#cd /Users/rshang/xmm_analysis/3HWC_J1928_p178/ID0902120101

cp /Users/rshang/xmm_analysis/analysis_code/*.sh .
cp /Users/rshang/xmm_analysis/analysis_code/*.py .
#
#
#
#
#
#sh extract_data.sh
#sh run_cifbuild.sh
#sh clean.sh 
#sh run_ana_chains.sh

#sh run_rename.sh

sh run_filtering.sh
#
#sh run_src_detect.sh

#sh run_create_regions.sh
#python3 convert_regions.py
#
#sh run_extract_spec_from_region.sh

#sh run_specgroup.sh
#sh run_img_prod.sh

