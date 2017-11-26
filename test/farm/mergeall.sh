outputDir=/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8

for fileFull in $outputDir/*
do
    ../merge $fileFull/convMats
    #echo $fileFull/convMats
    mv $fileFull/convMats/conv_F2.h5  $fileFull/conv_F2.h5
    mv $fileFull/convMats/conv_FL.h5  $fileFull/conv_FL.h5

done
