#! /bin/csh

set inputDir = /star/data19/embedding/pp500_production_2011/Psi2sEE_$1_20142601
#set inputDir = /eliza17/star/starprod/embedding/production2011_500GeV/$1_$2_2012815
if (-e datalist_$1.list) rm datalist_$1.list 
touch datalist_$1.list
echo $inputDir

foreach f (`find $inputDir/*/ -name 'st_physics_adc*'`)
 foreach put (`find $f -name '*.event.root'`)
  echo "$put">>datalist_$1.list
  set baseName = `basename $put`
  echo "$baseName"
end
end
