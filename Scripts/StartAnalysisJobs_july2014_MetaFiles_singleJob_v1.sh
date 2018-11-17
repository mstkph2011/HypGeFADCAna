#!/bin/bash

### script to start jobs on himster for germanium psa using go4
### v1.1 not tested! 			changed walltime for high sigmaGaus

SecondRun=0

SubPath=july2014						### subpath of subdirectory for analysis output

user=$USER				#### this is taken from system variable and used for job sending via the "double queue" (wait with job submissing if the internal himster queue is to full and send the jobs when there is enough space)

MWDmin=200				###Optimum 200
MWDmax=200
MWDstep=10
#if only a fixed value for sigma gaus should be used make min = max
sigmaGausmin=3
sigmaGausmax=3
sigmaGausstep=2
#if only a fixed value for sigma bil should be used make min = max
sigmaBilmin=0
#100
sigmaBilmax=0
#900
sigmaBilstep=20

AnaLibDir=/home/${user}/work/HypGeFADCAna								### path to analysis library


echo $DataInputFilePath
#parameters of GO4 analysis
### MWDm taken from loop values, see above for values

#MAl=120
MALmin=100
MALmax=100
MALstep=10

NumberOfSmoothings=100				### only used for rectangular or weighted average filter
FilterWidth=3
### sigmaGaus and Bil taken from loop values, see above for values
#tau=5383;							####z.Z. zwischen 6200 und 6230
taumin=6210;						###6210 aus pol2 fit
taumax=6210;	
taustep=100;	

EnableMA=1		
FilterType=0									### 0 = none, 1 = rectanglur, 2 = weighted average, 3 = gausian filter, 4 = bil filter
EnableBaselineCorrection=1

InputFile=/data/work/kpha1/steinen/COSYBeamtest/dataJuly2014/AnalysisInputMetaFiles/MetaFilesList.txt
echo blabla
while read line
#for ((Date=${Date}; Date<=${DateMax}; Date=$(($Date+${DateStep}))))
do
#echo $line
	InputFileAdd=$(echo $line | awk -F'/' '{print $9}')
	InputFileAdd=$(echo $InputFileAdd | awk -F'.' '{print $1}')
	#echo $InputFileAdd
	DataSubPath=dataJuly2014/${InputFileAdd}					### subpath of input data
	DataInputFilePath=${COSYTESTDATADIR}/${DataSubPath}							### input data directory

 SubDir=${COSYTESTANADIR}/${SubPath}				### complete path of subdirectory for analysis output
 if [ ! -d $SubDir ]; then 
   mkdir -p $SubDir
 fi

 JobLogPath=${COSYTESTANADIR}/${SubPath}/joblogs
 if [ ! -d $JobLogPath ]; then 
   mkdir -p $JobLogPath
 fi

 SimLogPath=${COSYTESTANADIR}/${SubPath}/analogs
 if [ ! -d $SimLogPath ]; then 
   mkdir -p $SimLogPath
 fi

 JobPath=${COSYTESTANADIR}/${SubPath}/jobs
 if [ ! -d $JobPath ]; then 
   mkdir -p $JobPath
 fi
 RunPath=${COSYTESTANADIR}/${SubPath}/runs
 if [ ! -d $RunPath ]; then 
   mkdir -p $RunPath
 fi
 jobcounter=1



	for ((MWDm=${MWDmin}; MWDm<=${MWDmax}; MWDm=$(($MWDm+${MWDstep}))))
	do
		for ((MAl=${MALmin}; MAl<=${MALmax} ; MAl=$(($MAl+${MALstep}))))
		do
			for ((tau=${taumin}; tau<=${taumax}; tau=$(($tau+${taustep}))))
			do
				for ((sigmaGaus=${sigmaGausmin}; sigmaGaus<=${sigmaGausmax}; sigmaGaus=$(($sigmaGaus+${sigmaGausstep}))))
				do
					for ((sigmaBil=${sigmaBilmin}; sigmaBil<=${sigmaBilmax}; sigmaBil=$(($sigmaBil+${sigmaBilstep}))))
					do
						
						SubSubDir=COSY_Ana_${InputFileAdd}___${MWDm},${MAl},${FilterType},${sigmaGaus},${sigmaBil},${tau}

						mkdir -p ${SubDir}/${SubSubDir}
		
						fileAdd=_${SubSubDir}
						ParameterFile=${COSYTESTANADIR}/july2014/DatabaseFirstAnalysisStep/ParametersFirstAnaStep${SubSubDir},MA.root
																																																							
						if [ ${SecondRun} -gt 0 ]; then
							fileAdd=${fileAdd}_SR
						fi		
						
if [ ! -f $JobPath/job_COSYJuly2014Ana_${fileAdd}.sh ]
then
				cat >$JobPath/job_COSYJuly2014Ana_${fileAdd}.sh <<EOF
#!/bin/bash

#PBS -j oe
#PBS -o ${JobLogPath}/job\${qfileAdd}.log
#PBS -V
### #PBS -l nodes=1:ppn=1,walltime=02:00:00,mem=4000mb
 #PBS -l nodes=1:ppn=1,walltime=02:00:00

export PATH=\${PBS_O_PATH}
cd \${PBS_O_WORKDIR}

echo go4analysis -number 0 -lib \${qRunPath}/run\${qfileAdd}/libGo4UserAnalysis.so -file @\${qline} -asf \${qSubDir}/\${qSubSubDir}/Ana\${qfileAdd}.root -x \${MWDm} \${MAl} \${NumberOfSmoothings} \${FilterWidth} \${sigmaGaus} \${sigmaBil} \${tau} \${EnableMA} \${FilterType} \${EnableBaselineCorrection} \${SecondRun} \${qParameterFile}

rm \${qSubDir}/\${qSubSubDir}/Ana\${qfileAdd}.root
go4analysis -number 0 -lib \${qRunPath}/run\${qfileAdd}/libGo4UserAnalysis.so -file @\${qline} -asf \${qSubDir}/\${qSubSubDir}/Ana\${qfileAdd}.root -x \${MWDm} \${MAl} \${NumberOfSmoothings} \${FilterWidth} \${sigmaGaus} \${sigmaBil} \${tau} \${EnableMA} \${FilterType} \${EnableBaselineCorrection} \${SecondRun} \${qParameterFile}
rm \${qRunPath}/run\${qfileAdd}/libGo4UserAnalysis.so

EOF
fi
						echo "jobcount: ${jobcounter}"
						jobcounter=$(($jobcounter+1))
	### submit job to batch system
						mkdir -p ${RunPath}/run$fileAdd
						cp ${AnaLibDir}/libGo4UserAnalysis.so ${RunPath}/run${fileAdd}/libGo4UserAnalysis.so
						#qsub -N job${fileAdd} -v qRunPath=\"${RunPath}\",qfileAdd=\"${fileAdd}\",qline=\"${line}\",qSubDir=\"$SubDir\",qSubSubDir=\"$SubSubDir\",qParameterFile=\"$ParameterFile\" $JobPath/job_COSYJuly2014Ana_${fileAdd}.sh
						#echo qsub -N job${fileAdd} -v qRunPath=\"${RunPath}\",qfileAdd=\"${fileAdd}\",qline=\"${line}\",qSubDir=\"$SubDir\",qSubSubDir=\"$SubSubDir\",qParameterFile=\"$ParameterFile\" $JobPath/job_COSYJuly2014Ana_${fileAdd}.sh
			go4analysis -number 0 -lib ${RunPath}/run${fileAdd}/libGo4UserAnalysis.so -file @${line} -asf ${SubDir}/${SubSubDir}/Ana${fileAdd}.root -x ${MWDm} ${MAl} ${NumberOfSmoothings} ${FilterWidth} ${sigmaGaus} ${sigmaBil} ${tau} ${EnableMA} ${FilterType} ${EnableBaselineCorrection} ${SecondRun} ${ParameterFile}

	### "double queue": check if there is enough space in the queue to send more jobs, if not wait 2 minutes and check again

						x=$(qstat | grep -c ${user})		
						while [ $x -gt 890 ]
						do
							echo "Still some jobs running: waiting to start new jobs"
							date
							sleep 2m #30s
							#x=$(ls | grep -c ext)
							x=$(qstat | grep -c ${user})
						done			### end of while "double queue" loop
					done 		### end of sigmaBil loop
				done		### end of sigmaGaus loop
			done		### end of tau loop
		done		###end of MAL loop
	done		### end of MWDm loop
done	< /data/work/kpha1/steinen/COSYBeamtest/dataJuly2014/AnalysisInputMetaFiles/MetaFilesList.txt #${InputFile}
#${go4Command}
