export OMP_NUM_THREADS=1

for fgtype in ETG_NE2001 ETG_YMW16 LTG_NE2001 LTG_YMW16 ALG_NE2001 ALG_YMW16
do
	echo "mpiexec.hydra -n 120 ./mcmc_frb.py -f frb_cat.txt -g ${fgtype} -o ${fgtype}"
    mpiexec.hydra -n 120 ./mcmc_frb.py -f frb_cat.txt -g ${fgtype} -o ${fgtype}
    mpiexec.hydra -n 120 ./mcmc_frb.py -f frb_cat.txt -g ${fgtype} -halo -o ${fgtype}_halo
	#mpiexec.hydra -n 60 ./mcmc_frb.py -f frb_cat.txt -g ${fgtype} -halo -upper -o ${fgtype}_upper_halo
	#mpiexec.hydra -n 60 ./mcmc_frb.py -f frb_cat.txt -g ${fgtype} -upper -o ${fgtype}_upper
done
