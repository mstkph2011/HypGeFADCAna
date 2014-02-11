#!/bin/bash


DataDir=${COSYTESTANADIR}/COSY
SubDir=${DataDir}/CombinedData
TxtDir=${COSYTESTANADIR}/COSY/txtfiles

rm ${TxtDir}/Filestofit.txt

ls ${SubDir}/*.root >> ${TxtDir}/Filestofit.txt
