
cd ODF/
rm *.TAR
rm *.FIT
rm *.ASC
rm *.SAS
rm MANIFEST.*
file=$(ls *.tar.gz)
#file=$(ls *.tar)
tar -xvf $file
file=$(ls *.TAR)
tar -xvf $file
exit 0
