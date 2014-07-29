#! /bin/csh

#positron
#set inputDir  = ./production2009_200GeV/$1_$2_20104802
#electron
set inputDir  = ./out_$1_$2/
set outputName = ./EffOutput/psiEff_2011_$1_$2
set logDir = ./logEff_$1_$2
mkdir $outputDir
mkdir $logDir
mkdir scriptEff

cp run.con EffrunAll_$2.job
echo $inputDir
cp ./run.csh ./run_tmp.csh
echo "root4star -b <<EOF">>run_tmp.csh
echo ".O2">>run_tmp.csh
echo -n '.x psiEff.C("'>>run_tmp.csh
echo -n $inputDir>>run_tmp.csh
echo -n '","'>>run_tmp.csh
#       echo -n $baseName>>run_tmp.csh
#       echo -n '","'>>run_tmp.csh
echo -n $outputName>>run_tmp.csh
#       echo -n '","'>>run_tmp.csh
echo -n '",'>>run_tmp.csh
echo -n $3>>run_tmp.csh
echo -n ','>>run_tmp.csh
echo -n $4>>run_tmp.csh
echo   ')'>>run_tmp.csh
echo ".q">>run_tmp.csh
echo "EOF">>run_tmp.csh
#
mv run_tmp.csh ./scriptEff/psiEff_run_$1_$2_$3_$4.csh
echo "Executable       = scriptEff/psiEff_run_$1_$2_$3_$4.csh">>EffrunAll_$2.job
echo "Output           = $logDir/psiEff_run_$1_$2_$3_$4.out">>EffrunAll_$2.job
echo "Error            = $logDir/psiEff_run_$1_$2_$3_$4.err">>EffrunAll_$2.job
echo "Log              = $logDir/psiEff_run_$1_$2_$3_$4.olog">>EffrunAll_$2.job
echo  "Queue" >>EffrunAll_$2.job
echo  "     " >>EffrunAll_$2.job

condor_submit EffrunAll_$2.job
