#!/bin/sh                                                                                                                       

gfixs=$HOME/gfix_from_gmat

iCdir=/home/LATTICE/fellow/kesasaki/confs/Set0/Kud013710Ks013710
GFdir=/home/LATTICE/fellow/kesasaki/confs/Set0/GF.dir
oCdir=$HOME/data/conf

for conf in `ls $GFdir`; do
    $gfixs 16 32 $iCdir/$conf $GFdir/$conf/GFmat.$conf $oCdir/$conf.gfix
done
