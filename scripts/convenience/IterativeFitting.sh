#!/bin/sh

if [ $# -lt 4 ]; then
    echo "`basename $0` [ifile] [ofile] [#.iter] [initial params]"
    exit 1
fi

######### ! Change Here ! #########
#Fit=$HOME/analysis_code-set/bin/fit
#Arg=$HOME/analysis_code-set/arg/arguments_fitting.txt
Fit=/Users/miiya/Dropbox/programs/analysis_code-set/bin/fit
Arg=/Users/miiya/Dropbox/programs/analysis_code-set/arg/arguments_fitting.txt

ifile=$1
ofile=$2

r_min=0
r_max=70
fit_func=4G

Arguments="-ifile $ifile -ofile $ofile -r_min $r_min -r_max $r_max -fit_func $fit_func"
###################################

Niter=$3
Args=($@)
params=""
for i in `seq 3 1 $(($#-1))`; do
    params="$params ${Args[$i]}"
done

echo "Initial parameters:$params"

for i in `seq 0 1 $Niter`; do
    TMP_A=(`$Fit -f $Arg $Arguments -ofile /dev/null -param $params @ | \
	grep -e "@ parameter infomations." -e "@ Fitting Results  : chisq/dof" -A 1 | \
	cut -d " " -f 3-`)
    chisq="${TMP_A[5]} ${TMP_A[6]} ${TMP_A[7]}"
    params=""
    for j in `seq 17 1 $((${#TMP_A[@]}-1))`; do
	params="$params ${TMP_A[$j]}"
    done
    echo "Iter=`printf %03d $((i+1))`: ChiSq = $chisq, Params =$params"
done

echo

$Fit -f $Arg $Arguments -param $params @
