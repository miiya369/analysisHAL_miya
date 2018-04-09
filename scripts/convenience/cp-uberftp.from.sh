#!/bin/sh

if [ $# -lt 2 ]; then
    echo "usage: `basename $0` [File(s)/Directory(s) (Path from gfarm/...)] [Dst. on local]"
    exit 1
fi

# For KEK
readonly CLIENT=scjldg01.sc.kek.jp

# For RIKEN
#readonly CLIENT=jldgriken01

readonly args=($@)
readonly Narg=$#

readonly CopyTo=${args[$((Narg-1))]}

ExitCode=0

############################
######### Function #########
function CpUberFTP() {
    local File=$1
    local Base=$2
    
    uberftp gsiftp://$CLIENT/$File file:$CopyTo/$Base
    if [ $? -eq 0 ]
    then
	echo "@   Successful to copy: $File"
    else
	echo "@   ERROR while copy  : $File, ExitCode = $ExitCode"
	ExitCode=1
    fi
}
######### Function #########
############################

echo "@@@ COPY FROM JLDG START: `date`"
echo

for i in `seq 0 1 $((Narg-2))`; do
    Fname=${args[$i]}
    BaseN=`dirname $Fname`
    
    lst_type=(`uberftp -ls -r gsiftp://$CLIENT/$Fname | awk '{print $1}'`)
    lst_file=(`uberftp -ls -r gsiftp://$CLIENT/$Fname | awk '{print $8}'`)
    
    if [ ${#lst_type[@]} -ne ${#lst_file[@]} ]; then
	echo "@@@ ERROR: Something happend. (1)"
	ExitCode=1
	continue
    fi
    
    for i in `seq 0 1 $((${#lst_type[@]}-1))`; do
	BASE=`dirname ${lst_file[$i]} | sed -e "s!$BaseN!!g"`
	if [ `echo ${lst_type[$i]} | grep d | wc -w` -eq 0 ]; then
	    CpUberFTP ${lst_file[$i]} $BASE
	elif [ `echo ${lst_type[$i]} | grep d | wc -w` -eq 1 ]; then
	    mkdir -p $CopyTo/$BASE/`basename ${lst_file[$i]}`
	else
	    echo "@@@ ERROR: Something happend. (2)"
            ExitCode=1
	fi
    done
done

echo
echo "@@@ COPY FROM JLDG  END : `date`"

exit $ExitCode
