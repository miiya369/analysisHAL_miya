#!/bin/sh

############# CHANGE HERE ##############
# Note: The name of variable should be 'x'

### 3G       (c.ens2; t=10; r=0.001-2.0fm; lim=1e-7; chisq/dof=0.811394)
#FitFunc="a1 * exp(-(x/b1)**2) + a2 * exp(-(x/b2)**2) + a3 * exp(-(x/b3)**2)"
#ParamsInit=(1524.09767404624 0.116916512674567 461.216801509866 0.299118278078935 -51.1603131391376 0.773100797986594)
#ParamsList=(a1 b1 a2 b2 a3 b3)

### 4G       (c.ens1; t=13; r=0.001-2.0fm; lim=1e-7; chisq/dof=)
#FitFunc="a1 * exp(-(x/b1)**2) + a2 * exp(-(x/b2)**2) + a3 * exp(-(x/b3)**2) + a4 * exp(-(x/b4)**2)"
#ParamsInit=(1524.09767404624 0.116916512674567 461.216801509866 0.299118278078935 -51.1603131391376 0.473100797986594 -11.1603131391376 0.773100797986594)
#ParamsList=(a1 b1 a2 b2 a3 b3 a4 b4)

### 2G1Ysq   (c.ens2; t=10; r=0.001-2.0fm; lim=1e-7; chisq/dof=0.426251)
#FitFunc="a1 * exp(-(x/a2)**2) + a3 * exp(-(x/a4)**2) + a5 * (1-exp(-a6*x**2))**2 * exp(-2*a7*x)/x**2"
#ParamsInit=(1276.655689454 0.100340021192985 860.055111519572 0.477599134464085 -25413.4438714745 1.23009311522288 3.1992278938778)
#ParamsList=(a1 a2 a3 a4 a5 a6 a7)

### 3G1Ysq   (c.ens2; t=10; r=0.001-2.0fm; lim=1e-7; chisq/dof=0.42635)
#FitFunc="a1 * exp(-(x/a2)**2) + a3 * exp(-(x/a4)**2) + a5 * exp(-(x/a6)**2) + a7 * (1-exp(-a8*x**2))**2 * exp(-2*a9*x)/x**2"
#ParamsInit=(1266.97646045365 0.100323252628672 51.6088989642378 0.304724368638692 819.581558388257 0.478065771766971 -43219.7200792367 0.964748841523314 3.34056801957467)
#ParamsList=(a1 a2 a3 a4 a5 a6 a7 a8 a9)

### 2G1Y     (c.ens2; t=10; r=0.001-2.0fm; lim=1e-7; chisq/dof=0.967228)
FitFunc="a1 * exp(-(x/a2)**2) + a3 * exp(-(x/a4)**2) + a5 * (1-exp(-a6*x**2)) * exp(-a7*x)/x"
ParamsInit=(1497.02243954189 0.118150020065909 444.329721588973 0.311776311948295 -588.246834072043 0.922070729303004 3.64571522933845)
ParamsList=(a1 a2 a3 a4 a5 a6 a7)

### 2G1Ysq1Y (c.ens2; t=10; r=0.001-2.0fm; lim=1e-7; chisq/dof=0.427729)
#FitFunc="a1 * exp(-(x/a2)**2) + a3 * exp(-(x/a4)**2) + a5 * (1-exp(-a6*x**2))**2 * exp(-2*a7*x)/x**2 + a8 * (1-exp(-a9*x**2)) * exp(-a10*x)/x"
#ParamsInit=(1275.40427582237 0.100346397539281 863.689182410647 0.477499374596227 -24791.0749366166 1.24445106087646 3.19251064617811 -152.336305847677 0.482523719483371 8.17863464907994)
#ParamsList=(a1 a2 a3 a4 a5 a6 a7 a8 a9 a10)

xDataColumn=1
yDataColumn=2
eDataColumn=3

FitLimit="1e-7"
########################################

if [ ${#ParamsList[@]} -ne ${#ParamsInit[@]} ]; then
    echo "#.params and #.initial value are different."; exit 1
fi

if [ $# -ne 1 -a $# -ne 3 ]; then
    echo "usage: `basename $0` [plot data file] ( [x-axis min] [x-axis max] )"; exit 1
fi

ifile=$1
Nparam=${#ParamsList[@]}
iparam="${ParamsList[0]}"
iparam_init="${ParamsList[0]}=${ParamsInit[0]}"

for i in `seq 1 1 $((Nparam-1))`; do
    iparam="$iparam, ${ParamsList[$i]}"
    iparam_init="$iparam_init; ${ParamsList[$i]}=${ParamsInit[$i]}"
done

if [ $# -eq 1 ]; then    
    gnuplot \
	-e "set fit logfile '/dev/null'" \
	-e "FIT_LIMIT = $FitLimit" \
	-e "$iparam_init" \
	-e "f(x) = $FitFunc" \
	-e "fit f(x) '$ifile' u $xDataColumn:$yDataColumn via $iparam" \
	-e "print ''" \
	-e "print 'Fit_Params_Results : ', $iparam" \
	-e "plot '$ifile' using $xDataColumn:$yDataColumn:$eDataColumn w yer, f(x); pause -1"
else
    x_min=$2
    x_max=$3
    gnuplot \
        -e "set fit logfile '/dev/null'" \
	-e "FIT_LIMIT = $FitLimit" \
	-e "$iparam_init" \
	-e "f(x) = $FitFunc" \
        -e "fit [x=$x_min:$x_max] f(x) '$ifile' u $xDataColumn:$yDataColumn:$eDataColumn via $iparam" \
	-e "print ''" \
	-e "print 'Fit_Params_Results : ', $iparam" \
	-e "plot '$ifile' using $xDataColumn:$yDataColumn:$eDataColumn w yer, f(x); pause -1"
fi
