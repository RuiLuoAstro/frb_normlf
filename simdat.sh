#-----------------------------------------------------------------------------------------

Nfrb=100
for galaxy_type in ETG LTG_NE2001 LTG_YMW16 ALG_NE2001 ALG_YMW16
do
	echo $galaxy_type

alpha=-1.0
logls=43.0

flux_thre=1.0 #Jy
specwidth=1000 #MHz
outputfile=./simudat/simdat_${galaxy_type}_${logls}_${alpha}.txt

#python simu.py -alpha $alpha -logls $logls -ns $Nfrb -thre $flux_thre -dnu $specwidth -type $galaxy_type -o $outputfile &
#echo "python simu.py -alpha $alpha -logls $logls -ns $Nfrb -thre $flux_thre -dnu $specwidth -type $galaxy_type -o $outputfile"
#-----------------------------------------------------------------------------------------

# simufrb_2.txt alpha=1.0 logls=44.0 flux_thre=1.0 #Jy specwidth=1000 #MHz
#galaxy_type=ETG outputfile=simdat_${galaxy_type}_${logls}_${alpha}.txt

echo "python simufrb.py -alpha $alpha -logls $logls -ns $Nfrb -thre $flux_thre -dnu $specwidth -type $galaxy_type -o $outputfile"
python simufrb.py -alpha $alpha -logls $logls -ns $Nfrb -thre $flux_thre -dnu $specwidth -type $galaxy_type -o $outputfile &



#-----------------------------------------------------------------------------------------

# simufrb_3.txt
alpha=1.0
logls=44.0
flux_thre=1.0 #Jy
specwidth=1000 #MHz
#galaxy_type=ETG
outputfile=./simudat/simdat_${galaxy_type}_${logls}_${alpha}.txt

echo "python simufrb.py -alpha $alpha -logls $logls -ns $Nfrb -thre $flux_thre -dnu $specwidth -type $galaxy_type -o $outputfile"
python simufrb.py -alpha $alpha -logls $logls -ns $Nfrb -thre $flux_thre -dnu $specwidth -type $galaxy_type -o $outputfile &

done
exit


