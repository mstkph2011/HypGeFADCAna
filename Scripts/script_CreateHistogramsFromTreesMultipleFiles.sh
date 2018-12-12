#!/bin/bash

#./SpecAna/CreateFilelistToFit.sh

SubSubDir=COSYnewMogon
TxtDir=${COSYTESTANADIR}/${SubSubDir}/txtfiles
OutDir=$COSYTESTANADIR/${SubSubDir}/Fit/

i=0
#more ${TxtDir}/Filestofit.txt
while read line ;
do
	#echo $line
	OutFilePath=$(dirname "${line}")
	#OutFileName=${line##*/}
	OutFileName=$(basename "${line}")
	FitFile=${OutDir}/Fit${OutFileName#"Tree"}
	#echo ${OutFilePath}/Histos${OutFileName}
	OutFile=${OutFilePath}/Histos${OutFileName}
	
	# root -l -b -q HypGeAnaCreateHistogramsFromTree.C\(\"${line}\"\,\"${OutFile}\"\,\"\"\)
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 
	#root -l -q -b DoTreeBasedT10DeriMaxCorrection.C\(\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\)
	#  #first round of corrections available
	#root -l -b -q HypGeAnaCreateHistogramsFromTree.C\(\"${line}\"\,\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\)
	
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_T10DeriMaxEnergyNorm_2
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_T1090EnergyNorm_2
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_TDeriMax90EnergyNorm_2
#	root -l -q -b DoTreeBasedT10DeriMaxCorrection.C\(\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\)
	###second roud of corrections available
#	root -l -b -q HypGeAnaCreateHistogramsFromTree.C\(\"${line}\"\,\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\)
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_T10DeriMaxEnergyNorm_2 Corr2_T1090EnergyNorm_2
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_T10DeriMaxEnergyNorm_2 Corr2_TDeriMax90EnergyNorm_2
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_T1090EnergyNorm_2 Corr2_T10DeriMaxEnergyNorm_2
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_T1090EnergyNorm_2 Corr2_TDeriMax90EnergyNorm_2
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_TDeriMax90EnergyNorm_2 Corr2_T10DeriMaxEnergyNorm_2
#	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_TDeriMax90EnergyNorm_2 Corr2_T1090EnergyNorm_2
	
	
	#if [ "$i" -eq "6" ]
	#then
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_T10DeriMaxEnergyNorm_2 Corr2_TDeriMax90EnergyNorm_2 cut
	#./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 
	##./SpecAna/RealSpecAna ${OutFile} ${FitFile} jülich2 Corr1_TDeriMax90EnergyNorm_2 Corr2_T1090EnergyNorm_2
	#fi
	#exit 0;
((i++))
done < ${TxtDir}/TreeFiles.txt


