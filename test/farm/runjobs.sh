address=/afs/desy.de/user/z/zlebcr/h1/TMD/Krakow/bkEvol/test/

outputDir=/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8

for fileFull in $outputDir/*
do
    mkdir $fileFull/convMats
    file=`basename $fileFull`

for i in `seq -w 0 45`
do
name=bkEvol_${file}_${i}
#name=disc_${i}
qsub_file=$address/farm/sub/${name}.sub


echo '#!/bin/sh' > $qsub_file
echo "#$ -N ${name}" >>  $qsub_file
echo "#$ -o $address/farm/out/${name}.out" >> $qsub_file
echo "#$ -e $address/farm/err/${name}.err" >> $qsub_file

echo "cd \$TMPDIR" >> $qsub_file
echo "time $address/convol $i < $fileFull/config_${file}.ini " >> $qsub_file
#echo "$address/disc/a.out" >> $qsub_file


#echo "cp conv_??_*.h5   /nfs/dust/cms/user/zlebcr/Krakow/convMat/"  >> $qsub_file
echo "cp conv_??_*.h5   $fileFull/convMats"  >> $qsub_file

done

done
