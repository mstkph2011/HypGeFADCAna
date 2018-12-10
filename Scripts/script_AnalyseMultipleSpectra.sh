#!/bin/bash



while getopts ":c:f:" arg; do
  case $arg in
    c) corr1=$OPTARG;;
    f) corr2=$OPTARG;;
  esac
done



echo $corr1


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
	./SpecAna/RealSpecAna $line ${OutDir}/Fit${OutFileName} jülich2 ${corr1} ${corr2}
((i++))
done < ${TxtDir}/Filestofit.txt


