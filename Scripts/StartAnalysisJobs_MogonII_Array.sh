#!/bin/bash

### script to start jobs on himster for germanium psa using go4
### v1.1 not tested! 			changed walltime for high sigmaGaus

SubPath=COSYnewMogon

user=$USER				#### this is taken from system variable and used for job sending via the "double queue" (wait with job submissing if the internal himster queue is to full and send the jobs when there is enough space)

MWDmin=200
MWDmax=200
MWDstep=20
#if only a fixed value for sigma gaus should be used make min = max
sigmaGausmin=9
#sigmaGausmax=13
sigmaGausmax=9
sigmaGausstep=2
#if only a fixed value for sigma bil should be used make min = max
sigmaBilmin=2000
#100
sigmaBilmax=2000
#900
sigmaBilstep=20

AnaLibDir=/home/steinen/work/HypGeFADCAna
DataInputFilePath=${COSYTESTDATADIR}/data/3110/2
#NumberOfFiles=106
NumberOfFiles=1

#parameters of GO4 analysis
### MWDm taken from loop values, see above for values
MAl=100
NumberOfSmoothings=100				### only used for rectangular or weighted average filter
FilterWidth=3
### sigmaGaus and Bil taken from loop values, see above for values
tau=5383;	
EnableMA=0		
FilterType=3									### 0 = none, 1 = rectanglur, 2 = weighted average, 3 = gausian filter, 4 = bil filter
EnableBaselineCorrection=1

echo ${COSYTESTANADIR}

SubDir=${COSYTESTANADIR}/${SubPath}
echo $SubDir
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

echo folders created
for ((MWDm=${MWDmin}; MWDm<=${MWDmax}; MWDm=$(($MWDm+${MWDstep}))))
do
	for ((sigmaGaus=${sigmaGausmin}; sigmaGaus<=${sigmaGausmax}; sigmaGaus=$(($sigmaGaus+${sigmaGausstep}))))
	do
		for ((sigmaBil=${sigmaBilmin}; sigmaBil<=${sigmaBilmax}; sigmaBil=$(($sigmaBil+${sigmaBilstep}))))
		do
			if [ $FilterType -eq 3 ];
			then
				SubSubDir=COSY_Ana${MWDm},${FilterType},${sigmaGaus}
			else
				if [ $FilterType -eq 4 ];
				then
					SubSubDir=COSY_Ana${MWDm},${FilterType},${sigmaGaus},${sigmaBil}
				fi
			fi
			mkdir -p ${SubDir}/${SubSubDir}
			for ((nFile=1; nFile<=${NumberOfFiles}; nFile=$(($nFile+1))))
			do
				
				fileAdd=_${SubSubDir}_file${nFile}
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
#SBATCH  -J job${fileAdd}

#SBATCH -o ${JobLogPath}/job${fileAdd}.log
#SBATCH -export=all
#SBATCH -t 120
#SBATCH -A m2_him_exp 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --mem-per-cpu=2700
#####   #SBATCH --partition=himster2_exp 
#SBATCH --partition=devel

cd \$SLURM_SUBMIT_DIR

#echo "Start Analysis of COSY data"
go4analysis -lib ${RunPath}/run${fileAdd}/libGo4UserAnalysis.so -file $DataInputFile -asf ${SubDir}/${SubSubDir}/Ana${fileAdd}.root -x ${MWDm} ${MAl} ${NumberOfSmoothings} ${FilterWidth} ${sigmaGaus} ${sigmaBil} ${tau} ${EnableMA} ${FilterType} ${EnableBaselineCorrection} 



EOF
				

				echo "jobcount: ${jobcounter}"
				jobcounter=$(($jobcounter+1))
### submit job to batch system
				mkdir -p ${RunPath}/run$fileAdd
				cp ${AnaLibDir}/libGo4UserAnalysis.so ${RunPath}/run${fileAdd}/libGo4UserAnalysis.so
				sbatch $JobPath/job${fileAdd}.sh | tee Jobinfo.txt
			done					### end of loop over files

### "double queue": check if there is enough space in the queue to send more jobs, if not wait 2 minutes and check again

			#x=$(qstat | grep -c ${user})		
			#while [ $x -gt 890 ]
			#do
				#echo "Still some jobs running: waiting to start new jobs"
				#date
				#sleep 2m #30s
				##x=$(ls | grep -c ext)
				#x=$(qstat | grep -c ${user})
			#done			### end of while "double queue" loop
		done 		### end of sigmaBil loop
	done		### end of sigmaGaus loop
done		### end of MWDm loop



