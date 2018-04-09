#!/bin/sh

if test $# -ne 3; then
    echo "usage : sh $0 [Tshift] [Ens] [Nconf/job]"
    exit 1
fi

Tshift=$1
Ens=$2
Nconf_job=$3

if [ $Ens == "1" ]
then
    Nconf=399
#    jobc=("a" "b" "c")
    jobc=("e" "d" "f")
    Njobc=3
#    jobc=("a" "b" "c" "d" "e")
#    jobc=("f" "g" "h" "i" "j")
#    Njobc=5
elif [ $Ens == "2" ]
then
    Nconf=400
    jobc=("a" "b" "c" "d")
    Njobc=4
elif [ $Ens == "3" ]
then
    Nconf=450
#    jobc=("a" "b" "c" "d" "e")
    jobc=("f" "g" "h" "i" "j")
    Njobc=5
#    jobc=("e" "f" "g" "h" "i" "j")
#    Njobc=6
else
    exit 1
fi

TSHIFT=`printf %03d $Tshift`

printJS=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/PrintJobscript_ForPACS-CS.sh
odir=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/pacs-cs_ens$Ens.RHQ/Jobscript/Tshift$TSHIFT
makeGL=/bgdata4/scnfqcd/miyamoto/work_Bridge1.3+HF/script/MakeGL_PACS-CS.gfix.py
gfile_dir=/bghome/scnfqcd/tkymiya/gfile_list/pacs-cs_ens${Ens}_gfix

NJS=$(($Nconf/$Nconf_job))
RestJob=$(($Nconf%$Nconf_job))

if [ $RestJob -eq 0 ]
then
    tmp=1
else
    tmp=0
fi

for i in `seq 0 1 $((NJS-$tmp))`
do
    conf_s=$(($i*$Nconf_job+1))
    conf_e=$((($i+1)*$Nconf_job))
    
    if [ $conf_e -gt $Nconf ]
    then
	conf_e=$Nconf
    fi

    conf=`printf %03d $conf_s`-`printf %03d $conf_e`
    if [ ! -e $gfile_dir/gfile_list.pacs-cs_ens${Ens}_gfix_$conf ]
    then
	$makeGL $conf_s $conf_e $Ens > $gfile_dir/gfile_list.pacs-cs_ens${Ens}_gfix_$conf
    fi
    
    if [ ! -e $odir ]
    then
        mkdir -p $odir
    fi
    
    $printJS $Tshift $i ${jobc[$((i%Njobc))]} $conf_s $conf_e $Ens > $odir/Tshift${TSHIFT}_$i.job
done
