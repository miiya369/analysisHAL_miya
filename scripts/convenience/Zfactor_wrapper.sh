#!/bin/sh

function calc() {
    awk "BEGIN {print $1}"
    return $?
}
ExtZfac=$HOME/analysis_code-set/scripts/convenience/ExtractZfactor.py

###################
### Change Here ###
base=$HOME/work_Bridge1.3+HF/working.dir/pacs-cs/full/c_quark/ens3/results.bin045/Analysis/mas

Had1PS=$base/Lambda_CG05_CG05_PS_correlator_mom0_fit
Had1SS=$base/Lambda_CG05_CG05_SS_correlator_mom0_fit

Had2PS=$base/Sigma_CG05_CG05_PS_correlator_mom0_fit
Had2SS=$base/Sigma_CG05_CG05_SS_correlator_mom0_fit

mass1=1.118781
mass2=1.183502

ParamHad1PS=(`$ExtZfac $Had1PS $mass1 15 20`)
ParamHad1SS=(`$ExtZfac $Had1SS $mass1 10 13`)
ParamHad2PS=(`$ExtZfac $Had2PS $mass2 15 20`)
ParamHad2SS=(`$ExtZfac $Had2SS $mass2 10 15`)
###################
###################

Nparam=${#ParamHad1PS[@]}
if [ $Nparam -ne ${#ParamHad1SS[@]} -o $Nparam -ne ${#ParamHad2PS[@]} -o $Nparam -ne ${#ParamHad2SS[@]} ]; then
    echo "#.param is different."
    exit 1
fi

sum=0
sum2=0
for i in `seq 0 1 $((Nparam-1))`; do
    Zfac1=`calc "sqrt(2.0*$mass1) * ${ParamHad1PS[$i]} / sqrt(${ParamHad1SS[$i]})"`
    Zfac2=`calc "sqrt(2.0*$mass2) * ${ParamHad2PS[$i]} / sqrt(${ParamHad2SS[$i]})"`
    
    tmp_sum=`calc "$Zfac1 / $Zfac2"`
    
    sum=`calc "$sum + $tmp_sum"`
    sum2=`calc "$sum2 + ($tmp_sum * $tmp_sum)"`
done

sum=`calc "$sum / $Nparam"`
sum2=`calc "$sum2 / $Nparam"`

err=`calc "sqrt( ( ($Nparam-1) * ($sum2 - ($sum * $sum)) )**2 )"`

echo "$sum +/- $err"
