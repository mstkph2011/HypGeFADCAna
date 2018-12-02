#!/bin/bash

SubSubSubDir=June2014
#SubSubSubDir=July2014
SubSubDir=COSYnewMogon
DataDir=${COSYTESTANADIR}/${SubSubDir}/${SubSubSubDir}
SubDir=${DataDir}
#SubDir=${DataDir}
TxtDir=${COSYTESTANADIR}/${SubSubDir}/txtfiles

rm ${TxtDir}/Filestofit.txt

ls ${SubDir}/*.root >> ${TxtDir}/Filestofit.txt
ls ${SubDir}/*.root
echo ${TxtDir}/Filestofit.txt
