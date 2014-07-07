#!/bin/bash

SubDir=june2014
DataDir=${COSYTESTANADIR}/${SubDir}
SubDir=${DataDir}/CombinedData
TxtDir=${COSYTESTANADIR}/${SubDir}/txtfiles

rm ${TxtDir}/Filestofit.txt

ls ${SubDir}/*.root >> ${TxtDir}/Filestofit.txt
