#!/bin/bash

### script to find jobs that didn't produce any output file, new jobs are sent automatically

#DataDir=/data/work/kpha1/rittgen/analysis/COSY
SubDir=june2014
DataDir=${COSYTESTANADIR}/${SubDir}
DataSubDir=${DataDir}/CombinedData


TxtDir=${COSYTESTANADIR}/${SubDir}/txtfiles
rm -f ${TxtDir}/WrongFilesInFolder.txt														# removes log file for erroneous runs , uncomment to automate
mkdir -p ${DataDir}/CombinedData
ls -d ${DataDir}/COSY_Ana*/ &> ${TxtDir}/AnaFolderList.txt				# write all folders with runs (parameter configs) to file
errors=0
rm ${TxtDir}MissingFiles.txt
while read line
do
	line=${line%/}																									### remove '/' from end of line  (shortest pattern matching, %% for longest)
	line=${line##*/}																								### remove longest pattern matching "*/" from the front of line		(# = shortest, ## = longest)
	#echo $line	
	#echo ls ${DataDir}/${line} | wc -w														
	x=$(ls ${DataDir}/${line} | wc -w )														# get number of files in run
	if [ $x -ne 106 ]																								# check if right number of files, otherwise run name to file 
	then
		errors=1
		echo "Missing files in $line   $x"
		for ((j = 1; j<107; j++))
		do
			ls ${DataDir}/${line}/Ana_${line}_file${j}.root &> /dev/null
			if [ $? -eq 2 ]							### check exit status of ls in the line above, 2 = file not found				, if (file not found) -> write jobname into file
			then
				echo -e "\t\tjob_${line}_file${j}.sh is missing"
				echo job_${line}_file${j}.sh >> ${TxtDir}MissingFiles.txt
			fi 
		done
		echo $line >> ${TxtDir}/WrongFilesInFolder.txt								# file must be empty before!!!!  			 ### this shows the runs with errors
	#else
		#ls ${DataDir}/${line}/*.root | xargs hadd -f1 -T ${DataSubDir}/${line}.root				# lists all root files in the run folder and combine them into one file, xargs appends the output of the ls to the command that follows
		#echo "File ${DataSubDir}/${line}.root created"
	fi
done < ${TxtDir}/AnaFolderList.txt

if [ $errors -ne 0 ]
then
	echo "there were some errors, check \"WrongFilesInFolder.txt\" to see in which run"
else
	echo "No errors found :)"
fi

### start jobs without output file

./StartMissingjobs.sh
