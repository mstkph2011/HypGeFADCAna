#!/bin/bash

SubSubDir=june2014
DataDir=${COSYTESTANADIR}/${SubSubDir}
SubDir=${DataDir}/CombinedData
TxtDir=${COSYTESTANADIR}/${SubSubDir}/txtfiles

rm ${TxtDir}/Filestofit.txt

ls ${SubDir}/*.root >> ${TxtDir}/Filestofit.txt
