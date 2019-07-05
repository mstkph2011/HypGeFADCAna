#!/bin/bash

#./SpecAna/CreateFilelistToFit.sh

SubSubDir=COSYnewMogon
TxtDir=${COSYTESTANADIR}/${SubSubDir}/txtfiles
OutDir=$COSYTESTANADIR/${SubSubDir}/Fit/

PeakNumber=0
#PeakNumber=2

#month=June
month=JuneBefore
#month=July

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
	
	#SpectrumModel=jülich2
	SpectrumModel=co60
	if [ "$month" == "July" ]
	then
		SpectrumModel=jülich3
		if [ $i -eq 0 ]
		then
			SpectrumModel=co60
		fi
	fi
	#((i++))
	echo $SpectrumModel
	#continue
	root -l -b -q HypGeAnaCreateHistogramsFromTree.C\(\"${line}\"\,\"${OutFile}\"\,\"\"\,${PeakNumber}\)
	
	echo ./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} 
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} 
	exit 0
	root -l -q -b DoTreeBasedT10DeriMaxCorrection.C\(\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\,${PeakNumber}\)
	#  #first round of corrections available
	root -l -b -q HypGeAnaCreateHistogramsFromTree.C\(\"${line}\"\,\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\,${PeakNumber}\)
	
	
	
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T10DeriMaxEnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T1090EnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_TDeriMax90EnergyNorm_${PeakNumber}
	root -l -q -b DoTreeBasedT10DeriMaxCorrection.C\(\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\,${PeakNumber}\)
	###second roud of corrections available
	root -l -b -q HypGeAnaCreateHistogramsFromTree.C\(\"${line}\"\,\"${OutFilePath}/Histos${OutFileName}\"\,\"${FitFile}\"\,${PeakNumber}\)
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T10DeriMaxEnergyNorm_${PeakNumber} Corr2_T1090EnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T10DeriMaxEnergyNorm_${PeakNumber} Corr2_TDeriMax90EnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T1090EnergyNorm_${PeakNumber} Corr2_T10DeriMaxEnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T1090EnergyNorm_${PeakNumber} Corr2_TDeriMax90EnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_TDeriMax90EnergyNorm_${PeakNumber} Corr2_T10DeriMaxEnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_TDeriMax90EnergyNorm_${PeakNumber} Corr2_T1090EnergyNorm_${PeakNumber}
	./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T10DeriMaxEnergyNorm_${PeakNumber} Corr2_TDeriMax90EnergyNorm_${PeakNumber} cut
	
	#if [ "$i" -eq "6" ]
	#then
	#./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_T10DeriMaxEnergyNorm_${PeakNumber} Corr2_TDeriMax90EnergyNorm_${PeakNumber} cut
	#./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} 
	##./SpecAna/RealSpecAna ${OutFile} ${FitFile} ${SpectrumModel} Corr1_TDeriMax90EnergyNorm_${PeakNumber} Corr2_T1090EnergyNorm_${PeakNumber}
	#fi
	#exit 0;
((i++))
done < ${TxtDir}/TreeFiles${month}.txt


