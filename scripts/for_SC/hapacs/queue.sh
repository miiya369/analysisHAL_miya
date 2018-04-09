#!/bin/sh

debgR=`qstat|grep debugq|grep R|wc -l`
debgQ=`qstat|grep debugq|grep Q|wc -l`
debgH=`qstat|grep debugq|grep H|wc -l`
comqR=`qstat|grep comq|grep R|wc -l`
comqQ=`qstat|grep comq|grep Q|wc -l`
comqH=`qstat|grep comq|grep H|wc -l`
tcagR=`qstat|grep tcag|grep R|wc -l`
tcagQ=`qstat|grep tcag|grep Q|wc -l`
tcagH=`qstat|grep tcag|grep H|wc -l`

echo;
echo "          job queue status         "
echo " ----------------+-----+-----+-----"
echo "      class      |  R  |  Q  |  H  "
echo " ------+---------+-----+-----+-----"
echo "  BASE | debugq  | `printf %03d $debgR` | `printf %03d $debgQ` | `printf %03d $debgH`"
echo "       | comq    | `printf %03d $comqR` | `printf %03d $comqQ` | `printf %03d $comqH`"
echo " ------+---------+-----+-----+-----"
echo "  TCA  | tcag    | `printf %03d $tcagR` | `printf %03d $tcagQ` | `printf %03d $tcagH`"
echo " ------+---------+-----+-----+-----"
echo;
