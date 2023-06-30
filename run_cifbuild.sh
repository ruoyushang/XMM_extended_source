source set_env.sh
echo "SAS_CCF="
echo $SAS_CCF
echo "SAS_ODF="
echo $SAS_ODF

# Output file: ccf.cif

cifbuild withccfpath=no analysisdate=now category=XMMCCF calindexset=$SAS_CCF fullpath=yes
odfingest odfdir=$SAS_ODF outdir=$SAS_ODF

exit 0
