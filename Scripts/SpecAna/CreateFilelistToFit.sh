#!/bin/bash

#SubSubDir=june2014
SubSubDir=july2014
DataDir=${COSYTESTANADIR}/${SubSubDir}
SubDir=${DataDir}/CombinedData
#SubDir=${DataDir}
TxtDir=${COSYTESTANADIR}/${SubSubDir}/txtfiles

rm ${TxtDir}/Filestofit.txt

ls ${SubDir}/*.root >> ${TxtDir}/Filestofit.txt
ls ${SubDir}/*.root
echo ${TxtDir}/Filestofit.txt
