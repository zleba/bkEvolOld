address=/home/zlebcr/prog/bkEvol/automate
for i in 90 100 110 120 
do
   sed "s/currentAlphaS/0.$i/" $address/steers/config.ini > $address/steersRep/config_as$i.ini
   mkdir -p data/eq8/as$i
   cp $address/steersRep/config_as$i.ini  data/eq8/as$i
   #$address/../iter < $address/steersRep/config_as$i.ini
   scp -r data/eq8/as$i naf-cms:/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8
   #scp -r data/eq8/as$i/* naf-cms:/nfs/dust/cms/user/zlebcr/Krakow/storage/eq8/as$i
   rm -r data/eq8/as$i
done
