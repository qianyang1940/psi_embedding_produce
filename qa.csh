#! /bin/csh

#positron
#set inputDir  = ./production2009_200GeV/$1_$2_20104802
#electron
set inputDir  = ./out_$1_$2/
set outputName = ./qaOutput/psiqa_2011_$1_$2
set logDir = ./logQA_$1_$2
mkdir $outputDir
mkdir $logDir
mkdir scriptQA

cp run.con QArun_$2.job
echo $inputDir
cp ./run.csh ./run_tmp.csh
echo "root4star -b <<EOF">>run_tmp.csh
echo ".O2">>run_tmp.csh
echo -n '.x qa.C("'>>run_tmp.csh
echo -n $inputDir>>run_tmp.csh
echo -n '","'>>run_tmp.csh
echo -n $outputName>>run_tmp.csh
echo -n '",'>>run_tmp.csh
echo -n $3>>run_tmp.csh
echo -n ','>>run_tmp.csh
echo -n $4>>run_tmp.csh
echo   ')'>>run_tmp.csh
echo ".q">>run_tmp.csh
echo "EOF">>run_tmp.csh
#
mv run_tmp.csh ./scriptQA/QA_run_$1_$2_$3_$4.csh

echo "Executable       = scriptQA/QA_run_$1_$2_$3_$4.csh">>QArun_$2.job
echo "Output           = $logDir/QA_run_$1_$2_$3_$4.out">>QArun_$2.job
echo "Error            = $logDir/QA_run_$1_$2_$3_$4.err">>QArun_$2.job
echo "Log              = $logDir/QA_run_$1_$2_$3_$4.olog">>QArun_$2.job
echo  "Queue" >>QArun_$2.job
echo  "     " >>QArun_$2.job

condor_submit QArun_$2.job
#
#
