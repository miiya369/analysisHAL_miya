#!/bin/sh

if [ $# -lt 4 ]; then
    echo "usage: `basename $0` [ifname] [mass] [max base] [rho0] [RA] [aA] [range_a 1] [range_a 2] ..."
    exit 1
fi

###### Change Here ######
#SolveSch=/$HOME/analysis_code-set/scripts/python_code_set/main/main_Solve.Sch.GaussExp.py
SolveSch=/Users/miiya/Dropbox/programs/analysis_code-set/scripts/python_code-set/Main/main_Solve.Sch.GaussExp.py

lat_space=0.0907
#lat_space=0.0845808

Energy_level=0
max_r=0.01
#########################

args=($@)
ifile=${args[0]}
mass=${args[1]}
max_base=${args[2]}
rho0=${args[3]}
RA=${args[4]}
aA=${args[5]}

printf "# range_a = "
for i in `seq 6 1 $((${#args[@]}-1))`; do
    printf "%lf, " ${args[$i]}
done
echo

for Nbase in `seq 1 1 $max_base`; do
    printf "%d " $Nbase
    for i in `seq 6 1 $((${#args[@]}-1))`; do
	EigE=(`$SolveSch $ifile --mass $mass --lat_space $lat_space --max_r $max_r \
	    --range_a ${args[$i]} --Nbase $Nbase --Nret $((Energy_level+1)) 2>&1 \
	    --rho0 $rho0 --RA $RA --aA $aA \
	    | grep "Eigen energy\[$Energy_level\]" \
	    | cut -d " " -f 5-7`)
	
	if [ ${#EigE[@]} -eq 0 ]; then EigE=(0.0 0.0 0.0); fi
	
	printf "%e %e " ${EigE[0]} ${EigE[2]}
    done
    echo
done
