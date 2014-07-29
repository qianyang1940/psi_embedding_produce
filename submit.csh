#! /bin/csh

#positron
#set inputDir  = ./production2009_200GeV/$1_$2_20104802
#electron
set inputDir = /star/data19/embedding/pp500_production_2011/Psi2sEE_$2_20142601
#set inputDir = /eliza17/star/starprod/embedding/production2011_500GeV/$1_$2_2012815
set outputDir = out_$1_$2
set logDir = log_$1_$2
mkdir $outputDir
mkdir $logDir
mkdir script_${1}_${2}

echo $inputDir

set maxJobs = $3

@ i=0
@ p=0
@ s=0

cp run.con runAll_$2.job
foreach f (`find $inputDir/*/ -name 'st_physics_adc*'`)
@ i++
 foreach put (`find $f -name '*.event.root'`)
 @ p++
  echo " *** $i $p $s $maxJobs ***"
  echo $put
  set baseName = `basename $put`
  set log = `basename $put`
       
  if (-e $logDir/$baseName.out) rm $logDir/$baseName.out  
  if (-e $outputDir/$baseName) then
    echo "file already exist"
  else 	
    @ s++
		cp run.csh run_tmp.csh
       echo "root4star -b <<EOF">>run_tmp.csh
       echo ".O2">>run_tmp.csh
       echo -n '.x doEmcEmbedEvent.C(1e9,"'>>run_tmp.csh
       echo -n $put>>run_tmp.csh
       echo -n '","'>>run_tmp.csh
       echo -n $outputDir>>run_tmp.csh
       echo '",kFALSE)'>>run_tmp.csh
       echo ".q">>run_tmp.csh
       echo "EOF">>run_tmp.csh
       mv run_tmp.csh ./script_${1}_${2}/run_$baseName.csh

	   echo "Executable       = script_${1}_${2}/run_$baseName.csh">>runAll_$2.job
	   echo "Output           = $logDir/run_$baseName.out">>runAll_$2.job
	   echo "Error            = $logDir/run_$baseName.err">>runAll_$2.job
	   echo "Log              = $logDir/run_$baseName.olog">>runAll_$2.job
	   echo  "Queue" >>runAll_$2.job
	   echo  "     " >>runAll_$2.job

 end
 if ($s >= $maxJobs ) break;  
end 
	condor_submit runAll_$2.job
