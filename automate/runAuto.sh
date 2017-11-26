address=/home/zlebcr/prog/bkEvol/automate
remoteAddress=/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8
for i in 100 # 120 130 140
#for i in 120 130 140
do
    if [ 1 = 0 ] ; then

   rm -r data/eq8/*
   mkdir -p data/eq8
   cp $address/steersRep/config_as$i.ini  data/eq8/as$i
   #$address/../iter < $address/steersRep/config_as$i.ini
   scp   naf-cms:/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8/as$i/*.[hi]* data/eq8
   fi

   ../iter < $address/data/eq8/config_as${i}.ini
   #scp -r data/eq8/as$i/* naf-cms:/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8/as$i
done
