#!/bin/bash

#Beamtime=June2014
Beamtime=July2014



MWDm=200
MAl=100 
NumberOfSmoothings=1 
FilterWidth=1 
sigmaGaus=1
sigmaBil=1 
tau=5339
EnableMA=1
FilterType=0
EnableBaselineCorrection=1

SecondAnalysisRound=0
ParameterFileName=/data/work/kpha1/steinen/COSYBeamtestAna/june2014/DatabaseFirstAnalysisStep/ParametersFirstAnaStepCOSY_Ana_1306_run1_1_20___200,100,0,3,0,6210,MA.root
BaselineValue=1413.6
TreeFile=TreeFile.root



#paths
#$COSYTESTANADIR 
#COSYnewMogon/${Beamtime}
#COSYnewMogon/jobs
#COSYnewMogon/joblogs

RunPath=${COSYTESTANADIR}/COSYnewMogon/runs/
OutputPath=${COSYTESTANADIR}/COSYnewMogon/${Beamtime}
JobLogPath=${COSYTESTANADIR}/COSYnewMogon/joblogs
JobPath=${COSYTESTANADIR}/COSYnewMogon/jobs
AnaLibDir=~/work/HypGeFADCAna


i=0
while read line ;
do
if [ -z "$line$" ];
then 
echo emtpy string
fi
Outputname=COSY${Beamtime}Dataset${i}_${MWDm},${MAl},${FilterType},${tau}_SR${SecondAnalysisRound}

ParameterFileName=${COSYTESTANADIR}/COSYnewMogon/Fit/FitCOSY${Beamtime}Dataset${i}_${MWDm},${MAl},${FilterType},${tau}_SR$((SecondAnalysisRound-1)).root

TreeFile=TreeCOSY${Beamtime}Dataset${i}_${MWDm},${MAl},${FilterType},${tau}_SR${SecondAnalysisRound}

cat >$JobPath/job${Outputname}.sh <<EOF
#!/bin/bash
#
#SBATCH  -J D${i}${Beamtime}

#SBATCH -o ${JobLogPath}/job${Outputname}.log
#SBATCH -e ${JobLogPath}/job${Outputname}.log
#SBATCH --open-mode=append
#SBATCH -export=all
#SBATCH --time=08:00:00 
#SBATCH -A m2_him_exp 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.steinen@him.uni-mainz.de
### #SBATCH --mem-per-cpu=2700
#SBATCH --mem-per-cpu=1350
#SBATCH --partition=himster2_exp 
### #SBATCH --partition=devel

cd \$SLURM_SUBMIT_DIR

#echo "Start Analysis of COSY data"
go4analysis -lib ${RunPath}/run${Outputname}/libGo4UserAnalysis.so -file @${line} -asf ${OutputPath}/${Outputname}.root -x ${MWDm} ${MAl} ${NumberOfSmoothings} ${FilterWidth} ${sigmaGaus} ${sigmaBil} ${tau} ${EnableMA} ${FilterType} ${EnableBaselineCorrection} ${SecondAnalysisRound} ${ParameterFileName} ${BaselineValue} ${OutputPath}/${TreeFile}.root



EOF
echo $ParameterFileName
mkdir -p ${RunPath}/run${Outputname}
cp ${AnaLibDir}/libGo4UserAnalysis.so ${RunPath}/run${Outputname}/libGo4UserAnalysis.so
#if [ $i -eq 14 ] 
#then
echo sbatch $JobPath/job${Outputname}.sh
sbatch $JobPath/job${Outputname}.sh
#fi
#

#echo $line
#echo $Outputname

((i++))
done 	< $SIMDATADIR/COSYBeamtest/data${Beamtime}/AnalysisInputMetaFiles/MetaFilesList.txt
