#!/bin/sh

if [ $# -lt 2 ]; then
    echo "usage: `basename $0` [File(s)/Directory(s)] [Dst. on JLDG (Path from gfarm/...)]"
    exit 1
fi

# For KEK
readonly CLIENT=scjldg01.sc.kek.jp

# For RIKEN
#readonly CLIENT=jldgriken01

readonly args=($@)
readonly Narg=$#

readonly CopyTo=${args[$((Narg-1))]}

current=`pwd`
ExitCode=0

############################
######### Function #########
function CpUberFTP() {
    local File=$1
    
    uberftp file:$File gsiftp://$CLIENT/$CopyTo/$File
    if [ $? -eq 0 ]
    then
	echo "@   Successful to copy: $File"
    else
	echo "@   ERROR while copy  : $File, ExitCode = $ExitCode"
	ExitCode=1
    fi
}
function CpUberFTP_Dir_Recursive() {
    local Dname=$1
    
    uberftp -mkdir gsiftp://$CLIENT/$CopyTo/$Dname
    
    for fname in `ls $Dname`; do
	local FNAME=`basename $fname`
	local Fname=$Dname/$FNAME
	if [ -f $Fname ]; then
            CpUberFTP $Fname
	elif [ -d $Fname ]; then
	    CpUberFTP_Dir_Recursive $Fname
	else
	    echo "Not-regular file or directory: $Fname"
	    ExitCode=1
	fi
    done
}
######### Function #########
############################

echo "@@@ COPY TO JLDG START: `date`"
echo

for i in `seq 0 1 $((Narg-2))`; do
    cd $current
    cd `dirname ${args[$i]}`
    Fname=`basename ${args[$i]}`
    
    if [ ! -e $Fname ]; then
	echo "No such a file or directory: `dirname ${args[$i]}`/$Fname"
	ExitCode=1
	continue
    fi
    
    if [ -f $Fname ]; then
	CpUberFTP $Fname
    elif [ -d $Fname ]; then
	CpUberFTP_Dir_Recursive $Fname
    else
	echo "Not-regular file or directory: `dirname ${args[$i]}`/$Fname"
	ExitCode=1
    fi
done

echo
echo "@@@ COPY TO JLDG  END : `date`"

exit $ExitCode
