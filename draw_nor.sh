export OMP_NUM_THREADS=1

for fgtype in ETG_NE2001 ETG_YMW16 LTG_NE2001 LTG_YMW16 ALG_NE2001 ALG_YMW16
#for fgtype in ETG_NE2001_halo ETG_YMW16_halo LTG_NE2001_halo LTG_YMW16_halo ALG_NE2001_halo ALG_YMW16_halo

do
    echo "New plot of ${fgtype} has been finished."
    ./pltnorm.py -f ./mn_out/${fgtype} -o ./plots/normal/${fgtype}.eps
done
