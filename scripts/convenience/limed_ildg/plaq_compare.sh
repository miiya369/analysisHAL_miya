#!/bin/sh

idir_org=/tmp/miyamoto.jldg/disk.jldg/gfarm/public/ILDG/JLDG/PACS-CS/RCNF2+1/RC32x64_B1900/Kud01370000Ks01364000
idir_fix=/tmp/miyamoto.jldg/disk.jldg/gfarm/pacscs/BBILQCD/miyamoto/conf.gfix/pacs-cs/RC32x64_B1900/kappa_013700.013640

Nconf=`ls $idir_fix | wc -l`
conf=`ls $idir_fix`

echo "#.conf = $Nconf"
echo

for CONF in $conf
do
    echo "conf = $CONF"
    ./plaq_compare 32 64 $idir_org/$CONF $idir_fix/$CONF
done
