
dir=$(pwd)
export SAS_CCFPATH=/Users/rshang/xmm_analysis/CCF
export SAS_CCF="$dir"/analysis/ccf.cif
file=$(ls ODF/*.SAS)
export SAS_ODF="$dir"/"$file"

export ELOW='0'
export EHIGH='12000'
