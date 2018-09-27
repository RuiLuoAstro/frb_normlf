#!/bin/sh

for galaxy_type in ETG LTG_NE2001 LTG_YMW16 ALG_NE2001 ALG_YMW16
do
	echo $galaxy_type
#galaxy_type=ETG 
alpha=-1.0
logls=43.0

./pltnorm.py -f ./mn_out/simdat_${galaxy_type}_${logls}_${alpha} -o ./plots/simu/simdat_${galaxy_type}_${alpha}.eps -title "Mock data"

./pltupper.py -f ./mn_out/simdat_${galaxy_type}_${logls}_${alpha}_upper -o ./plots/simu/simdat_${galaxy_type}_${alpha}_upper.eps -title "Mock data"

#-----------------------------------------------------------------------------------------

alpha=1.0
logls=44.0

./pltnorm.py -f ./mn_out/simdat_${galaxy_type}_${logls}_${alpha} -o ./plots/simu/simdat_${galaxy_type}_${alpha}.eps -title "Mock data"

./pltupper.py -f ./mn_out/simdat_${galaxy_type}_${logls}_${alpha}_upper -o ./plots/simu/simdat_${galaxy_type}_${alpha}_upper.eps -title "Mock data"

done
exit
