export OMP_NUM_THREADS=1

for fgtype in ETG_NE2001_upper ETG_YMW16_upper LTG_NE2001_upper LTG_YMW16_upper ALG_NE2001_upper ALG_YMW16_upper
#for fgtype in ETG_NE2001_upper_halo ETG_YMW16_upper_halo LTG_NE2001_upper_halo LTG_YMW16_upper_halo ALG_NE2001_upper_halo ALG_YMW16_upper_halo

do
    echo "New plot of ${fgtype} has been finished."
    ./pltupper.py -f ./mn_out/${fgtype} -o ./plots/upper/${fgtype}.eps
done
