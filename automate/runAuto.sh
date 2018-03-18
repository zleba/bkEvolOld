tag=eq8gen
address=/home/zlebcr/prog/bkEvol/automate
remoteAddress=/nfs/dust/cms/user/zlebcr/Krakow/storage/$tag

#names="as100fr10kt0   as110fr07kt2  as110fr10kt0  as110fr10kt2  as110fr10kt4  as110fr13kt2  as120fr10kt0  as90fr10kt0"
names="as110fr10kt2"

for name in $names
#for i in 100 # 120 130 140
#for i in 120 130 140
do
    if [ 1 = 1 ] ; then

   rm -r data/$tag/*
   mkdir -p data/$tag
   #cp $address/steersRep/config_as$i.ini  data/$tag/as$i
   #$address/../iter < $address/steersRep/config_as$i.ini
   scp   naf-cms:/nfs/dust/cms/user/zlebcr/Krakow/storage/$tag/${name}/*.[hi]* data/$tag
   fi

   Input=`ls $address/data/$tag/*.ini`
   ../iter < $Input
   cp data/eq8gen/fit2Parm.dat results/${name}.dat
   #scp -r data/eq8/as$i/* naf-cms:/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8/as$i
done
