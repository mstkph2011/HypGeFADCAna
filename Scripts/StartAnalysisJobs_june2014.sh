#!/bin/bash

### script to start jobs on himster for germanium psa using go4
### v1.1 not tested! 			changed walltime for high sigmaGaus

SubPath=june2014						### subpath of subdirectory for analysis output
DataSubPath=dataJune2014/1006/run1					### subpath of input data

user=$USER				#### this is taken from system variable and used for job sending via the "double queue" (wait with job submissing if the internal himster queue is to full and send the jobs when there is enough space)

MWDmin=200				###Optimum 200
MWDmax=200
MWDstep=2
#if only a fixed value for sigma gaus should be used make min = max
sigmaGausmin=3
sigmaGausmax=3
sigmaGausstep=2
#if only a fixed value for sigma bil should be used make min = max
sigmaBilmin=2000
#100
sigmaBilmax=2000
#900
sigmaBilstep=20

AnaLibDir=/home/${user}/work/HypGeFADCAna								### path to analysis library
DataInputFilePath=${COSYTESTDATADIR}/${DataSubPath}							### input data directory
NumberOfFiles=20																				### number of input files
echo $DataInputFilePath
#parameters of GO4 analysis
### MWDm taken from loop values, see above for values

#MAl=120
MALmin=100
MALmax=100
MALstep=2

NumberOfSmoothings=100				### only used for rectangular or weighted average filter
FilterWidth=3
### sigmaGaus and Bil taken from loop values, see above for values
#tau=5383;	
taumin=6200;						###over 5300
taumax=7000;	
taustep=30;	

EnableMA=1		
FilterType=3									### 0 = none, 1 = rectanglur, 2 = weighted average, 3 = gausian filter, 4 = bil filter
EnableBaselineCorrection=1



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

#MWDm=${MWDmin}
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
					if [ $FilterType -ne 4 ];
					then
						SubSubDir=COSY_Ana${MWDm},${MAl},${FilterType},${sigmaGaus},${tau}
					else
						if [ $FilterType -eq 4 ];
						then
							SubSubDir=COSY_Ana${MWDm},${MAl},${FilterType},${sigmaGaus},${sigmaBil},${tau}
						fi
					fi
					mkdir -p ${SubDir}/${SubSubDir}
					for ((nFile=1; nFile<=${NumberOfFiles}; nFile=$(($nFile+1))))
					do
				
						fileAdd=_${SubSubDir}_${nFile}    ##_file${nFile}
						echo $fileAdd
						if [ ${nFile} -lt 10 ]; then
						DataInputFile=${DataInputFilePath}/data000${nFile}
						else
							if [ ${nFile} -lt 100 -a ${nFile} -gt 9 ]; 
							then
								DataInputFile=${DataInputFilePath}/data00${nFile}
							else
								if  [ ${nFile} -gt 99 -a ${nFile} -lt 1000 ];
								then
									DataInputFile=${DataInputFilePath}/data0${nFile}
								else
									if  [ ${nFile} -gt 999 ]; 
										then
										DataInputFile=${DataInputFilePath}/data${nFile}
									fi
								fi
							fi
						fi			

					cat >$JobPath/job${fileAdd}.sh <<EOF
#!/bin/bash
#
#PBS -N job${fileAdd}
#PBS -j oe
#PBS -o ${JobLogPath}/job${fileAdd}.log
#PBS -V
#PBS -l nodes=1:ppn=1,walltime=02:00:00

cd \$PBS_O_WORKDIR

#echo "Start Analysis of COSY data"
go4analysis -lib ${RunPath}/run${fileAdd}/libGo4UserAnalysis.so -file $DataInputFile -asf ${SubDir}/${SubSubDir}/Ana${fileAdd}.root -x ${MWDm} ${MAl} ${NumberOfSmoothings} ${FilterWidth} ${sigmaGaus} ${sigmaBil} ${tau} ${EnableMA} ${FilterType} ${EnableBaselineCorrection} &> ${SimLogPath}/ana${fileAdd}.log


EOF
				
				echo "jobcount: ${jobcounter}"
				#echo "go4analysis -lib ${RunPath}/run${fileAdd}/libGo4UserAnalysis.so -file $DataInputFile -asf ${SubDir}/${SubSubDir}/Ana${fileAdd}.root -x ${MWDm} ${MAl} ${NumberOfSmoothings} ${FilterWidth} ${sigmaGaus} ${sigmaBil} ${tau} ${EnableMA} ${FilterType} ${EnableBaselineCorrection} &> ${SimLogPath}/ana${fileAdd}.log"
				jobcounter=$(($jobcounter+1))
### submit job to batch system
				mkdir -p ${RunPath}/run$fileAdd
				cp ${AnaLibDir}/libGo4UserAnalysis.so ${RunPath}/run${fileAdd}/libGo4UserAnalysis.so
				qsub $JobPath/job${fileAdd}.sh | tee Jobinfo.txt
			done					### end of loop over files

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

