#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage : sh $0 [Path to Gfiles] [#.conf/job] [#.job]"
    exit 1
fi

#echo "Safty Device is working."
#exit 0

######### ! CHANGE HERE ! #########
PrintJS=/bghome/scnfqcd/tkymiya/bgdata4/script/PrintJobscript.sh
PrintGL=/bghome/scnfqcd/tkymiya/bgdata4/script/PrintGfileList.sh

JSdir=/bghome/scnfqcd/tkymiya/tmp/for_cp-pacs/bridge-1.3.2-hal.v01/src2/run.cp-pacs/JobScript
GLdir=/bghome/scnfqcd/tkymiya/tmp/for_cp-pacs/bridge-1.3.2-hal.v01/src2/run.cp-pacs/gfile.lst

BaseNameJS=run
BaseNameGL=gfile_list.cp-pacs+JLQCD_ens3.gfix
###################################

PathToConf=$1
Nconf=`ls $PathToConf | wc -w`

Nconf_job=$2
Njob=$3

for i in `seq 0 1 $((Njob-1))`
do
    confN_s=$((i       * Nconf_job + 1))
    confN_e=$((confN_s + Nconf_job - 1))
    if [ $confN_e -gt $Nconf ]
    then
	confN_e=$Nconf
    fi
    
    $PrintJS $Nconf $Nconf_job $i > $JSdir/${BaseNameJS}_$i.job
    $PrintGL $PathToConf $Nconf_job $i > $GLdir/${BaseNameGL}_`printf %03d $confN_s`-`printf %03d $confN_e`
done
