#!/bin/bash

if [ $# -eq 0 ]; then
    echo "usage: `basename $0` [FitFunc] [params] ..."
    exit 1
fi

Nparam=$(($#-1))
p=($@)
FitFunc=${p[0]}

case $FitFunc in
    3G      ) if [ $Nparam -ne 6  ]; then echo "ERROR: #.param is differ, exit."; exit; fi
	printf "(${p[1]})*exp(-(x/(${p[2]}))**2)+"
	printf "(${p[3]})*exp(-(x/(${p[4]}))**2)+"
	printf "(${p[5]})*exp(-(x/(${p[6]}))**2)"; echo;;
    2G1Ysq  ) if [ $Nparam -ne 7  ]; then echo "ERROR: #.param is differ, exit."; exit; fi
	printf "(${p[1]})*exp(-(x/(${p[2]}))**2)+"
	printf "(${p[3]})*exp(-(x/(${p[4]}))**2)+"
	printf "(${p[5]})*(1-exp(-(${p[6]})*x**2))**2*exp(-2*(${p[7]})*x)/x**2"; echo;;
    3G1Ysq  ) if [ $Nparam -ne 9  ]; then echo "ERROR: #.param is differ, exit."; exit; fi
	printf "(${p[1]})*exp(-(x/(${p[2]}))**2)+"
	printf "(${p[3]})*exp(-(x/(${p[4]}))**2)+"
	printf "(${p[5]})*exp(-(x/(${p[6]}))**2)+"
	printf "(${p[7]})*(1-exp(-(${p[8]})*x**2))**2*exp(-2*(${p[9]})*x)/x**2"; echo;;
    2G1Y    ) if [ $Nparam -ne 7  ]; then echo "ERROR: #.param is differ, exit."; exit; fi
	printf "(${p[1]})*exp(-(x/(${p[2]}))**2)+"
	printf "(${p[3]})*exp(-(x/(${p[4]}))**2)+"
	printf "(${p[5]})*(1-exp(-(${p[6]})*x**2))*exp(-(${p[7]})*x)/x"; echo;;
    2G1Ysq1Y) if [ $Nparam -ne 10 ]; then echo "ERROR: #.param is differ, exit."; exit; fi
	printf "(${p[1]})*exp(-(x/(${p[2]}))**2)+"
	printf "(${p[3]})*exp(-(x/(${p[4]}))**2)+"
	printf "(${p[5]})*(1-exp(-(${p[6]})*x**2))**2*exp(-2*(${p[7]})*x)/x**2+"
	printf "(${p[8]})*(1-exp(-(${p[9]})*x**2))*exp(-(${p[10]})*x)/x"; echo;;
    *       ) echo "ERROR: FitFunc=$FitFunc has not implemented yet, exit."; exit;;
esac
