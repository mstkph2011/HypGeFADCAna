#!/bin/bash

#this script starts missing jobs, use another script (FindMissingFiles.sh) before to Create the input file

DataDir=${COSYTESTANADIR}/COSY

TxtDir=${COSYTESTANADIR}/COSY/txtfiles
InputFile=${TxtDir}MissingFiles.txt
while read line
do
	qsub ${DataDir}/jobs/$line
done < $InputFile
