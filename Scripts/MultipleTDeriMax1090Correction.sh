#!/bin/bash

Beamtime=June2014
MaxFile=14
SecondRun=1
echo $MaxFile
for i in {0..14}
do
	
	echo root -l -b -q DoTDeriMax1090Correction.C'("/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/'${Beamtime}'/COSY'${Beamtime}'Dataset'${i}'_200,100,0,5339_SR'${SecondRun}'.root","/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSY'${Beamtime}'Dataset'${i}'_200,100,0,5339_SR'${SecondRun}'.root")'
	root -l -b -q DoTDeriMax1090Correction.C'("/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/'${Beamtime}'/COSY'${Beamtime}'Dataset'${i}'_200,100,0,5339_SR'${SecondRun}'.root","/lustre/miifs05/scratch/him-specf/hyp/steinen/COSYBeamtestAna/COSYnewMogon/Fit/FitCOSY'${Beamtime}'Dataset'${i}'_200,100,0,5339_SR'${SecondRun}'.root")'
done


