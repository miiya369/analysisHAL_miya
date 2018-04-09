#!/bin/sh

tcagR=`qstat|grep tmiyamoto|grep tcag|grep R|wc -l`
tcagQ=`qstat|grep tmiyamoto|grep tcag|grep Q|wc -l`
tcagH=`qstat|grep tmiyamoto|grep tcag|grep H|wc -l`
comqR=`qstat|grep tmiyamoto|grep comq|grep R|wc -l`
comqQ=`qstat|grep tmiyamoto|grep comq|grep Q|wc -l`
comqH=`qstat|grep tmiyamoto|grep comq|grep H|wc -l`
debugR=`qstat|grep tmiyamoto|grep debugq|grep R|wc -l`
debugQ=`qstat|grep tmiyamoto|grep debugq|grep Q|wc -l`
debugH=`qstat|grep tmiyamoto|grep debugq|grep H|wc -l`

echo ""
echo "  class  |  R  |  Q  |  H"
echo " --------+-----+-----+-----"
echo "  debugq | `printf %03d $debugR` | `printf %03d $debugQ` | `printf %03d $debugH`"
echo "  comq   | `printf %03d $comqR` | `printf %03d $comqQ` | `printf %03d $comqH`"
echo " --------+-----+-----+-----"
echo "  tcag   | `printf %03d $tcagR` | `printf %03d $tcagQ` | `printf %03d $tcagH`"
echo " --------+-----+-----+-----"
echo ""
