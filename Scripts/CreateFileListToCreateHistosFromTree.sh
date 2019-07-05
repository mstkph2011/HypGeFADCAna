#!/bin/bash

#SubSubSubDir=June2014
SubSubSubDir=July2014
SubSubDir=COSYnewMogon
DataDir=${COSYTESTANADIR}/${SubSubDir}/${SubSubSubDir}
SubDir=${DataDir}
#SubDir=${DataDir}
TxtDir=${COSYTESTANADIR}/${SubSubDir}/txtfiles

rm ${TxtDir}/Filestofit.txt

ls ${SubDir}/Tree*.root >> ${TxtDir}/TreeFiles.txt
ls ${SubDir}/Tree*.root
echo ${TxtDir}/TreeFiles.txt
