for galaxy_type in ETG LTG_NE2001 LTG_YMW16 ALG_NE2001 ALG_YMW16
do
	echo $galaxy_type

alpha=-1.0
logls=43.0

flux_thre=1.0 #Jy
bandwidth=1000 #MHz
#galaxy_type=ETG
outputfile=simdat_${galaxy_type}_${logls}_${alpha}
inputfile=./simudat/simdat_${galaxy_type}_${logls}_${alpha}.txt

export OMP_NUM_THREADS=2
mpiexec.hydra -n 128 python mcmc_simu.py -f $inputfile -o $outputfile -g ${galaxy_type}
#-----------------------------------------------------------------------------------------

alpha=-1.0
logls=43.0

flux_thre=1.0 #Jy
bandwidth=1000 #MHz
#galaxy_type=ETG
outputfile=simdat_${galaxy_type}_${logls}_${alpha}_upper
inputfile=./simudat/simdat_${galaxy_type}_${logls}_${alpha}.txt

export OMP_NUM_THREADS=2
mpiexec.hydra -n 128 python mcmc_simu.py -upper -f $inputfile -o $outputfile -g ${galaxy_type}
#-----------------------------------------------------------------------------------------


alpha=1.0
logls=44.0

flux_thre=1.0 #Jy
bandwidth=1000 #MHz
#galaxy_type=ETG
outputfile=simdat_${galaxy_type}_${logls}_${alpha}
inputfile=./simudat/simdat_${galaxy_type}_${logls}_${alpha}.txt

export OMP_NUM_THREADS=2
mpiexec.hydra -n 128 python mcmc_simu.py -f $inputfile -o $outputfile -g ${galaxy_type}
#-----------------------------------------------------------------------------------------

alpha=1.0
logls=44.0

flux_thre=1.0 #Jy
bandwidth=1000 #MHz
#galaxy_type=ETG
outputfile=simdat_${galaxy_type}_${logls}_${alpha}_upper
inputfile=./simudat/simdat_${galaxy_type}_${logls}_${alpha}.txt

export OMP_NUM_THREADS=2
mpiexec.hydra -n 128 python mcmc_simu.py -upper -f $inputfile -o $outputfile -g ${galaxy_type}
done
exit


