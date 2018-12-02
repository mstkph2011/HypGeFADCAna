#!/bin/bash

#./SpecAna/CreateFilelistToFit.sh

SubSubDir=COSYnewMogon
TxtDir=${COSYTESTANADIR}/${SubSubDir}/txtfiles
OutDir=$COSYTESTANADIR/${SubSubDir}/Fit/

i=1
#more ${TxtDir}/Filestofit.txt
while read line ;
do
	#echo $line
	OutFileName=${line##*/}
	#echo ${OutDir}/Fit${OutFileName}
	./SpecAna/RealSpecAna $line ${OutDir}/Fit${OutFileName}
((i++))
done < ${TxtDir}/Filestofit.txt


