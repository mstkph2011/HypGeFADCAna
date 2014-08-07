#!/bin/bash

### script to combine the output of single runs to a combined root file
Sub=june2014
DataDir=${COSYTESTANADIR}/${Sub}
#DataDir=/data/work/kpha1/rittgen/analysis/COSYnoNoise

SubDir=${DataDir}/CombinedData
TxtDir=${COSYTESTANADIR}/${Sub}/txtfiles
NumberOfFiles=20

SubPath=june2014
JobLogPath=${COSYTESTANADIR}/${SubPath}/joblogs
if [ ! -d $JobLogPath ]; then 
  mkdir -p $JobLogPath
fi
SimLogPath=${COSYTESTANADIR}/${SubPath}/analogs
if [ ! -d $SimLogPath ]; then 
  mkdir -p $SimLogPath
fi
JobPath=${COSYTESTANADIR}/${SubPath}/jobs
if [ ! -d $JobPath ]; then 
  mkdir -p $JobPath
fi


rm -f ${DataDir}/WrongFilesInFolder.txt														# removes log file for erroneous runs , uncomment to automate
mkdir -p ${SubDir}
ls -d ${DataDir}/COSY_Ana*/ &> ${TxtDir}/AnaFolderList.txt				# write all folders with runs (parameter configs) to file
errors=0
while read line
do
	line=${line%/}																									# remove '/' from end of line
	line=${line##*/}
	echo $line	
	echo ls ${DataDir}/${line} | wc -w													
	x=$(ls ${DataDir}/${line} | wc -w )														# get number of files in run
	if [ $x -ne ${NumberOfFiles} ]																								# check if right number of files, otherwise run name to file 
	then
		errors=1
		echo $line
		echo $line >> ${TxtDir}/WrongFilesInFolder.txt								# file must be empty before!!!!
	else
		cat >${JobPath}/jobCombo${line}.sh <<EOF
		#!/bin/bash
		#
		#PBS -N job${line}
		#PBS -j oe
		#PBS -o ${JobLogPath}/job${line}.log
		#PBS -V
		#PBS -l nodes=1:ppn=1,walltime=02:00:00

		cd \$PBS_O_WORKDIR

		#echo "Combo job"
		ls ${DataDir}/${line}/*.root | xargs hadd -f5 -T ${SubDir}/${line}.root				# lists all root files in the run folder and combine them into one file, xargs appends the output of the ls to the command that follows

EOF
	
		#ls ${DataDir}/${line}/*.root | xargs hadd -f5 -T ${SubDir}/${line}.root				# lists all root files in the run folder and combine them into one file, xargs appends the output of the ls to the command that follows
		qsub ${JobPath}/jobCombo${line}.sh
		echo "File ${SubDir}/${line}.root created"
	fi
done < ${TxtDir}/AnaFolderList.txt

if [ $errors -ne 0 ]
then
	echo "there were some errors, check \"WrongFilesInFolder.txt\" to see in which run"
fi
