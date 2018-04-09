#!/bin/sh

if test $# -ne 2
then 
    echo " usage : sh $0 [i/h] [conf]"
    exit -1
fi

set -e

SET=$1
IDIR=/work/HALQCD/tmiyamoto/run-${SET}/bin/DownLoad/CONF/run-${SET}
ODIR=/work/HALQCD/tmiyamoto
CONF=`printf %06d $2`
PREC=1.e-14
ITER=10000000

KONV96=/work/HALQCD/tmiyamoto/working.dir/gfix.Kconf_ha-pacs/work/konv96
CoulombGF=/work/HALQCD/tmiyamoto/working.dir/gfix.Kconf_ha-pacs/work/CoulombGaugeFixingSU3_4D_96x13
SCRIPT=/home/HALQCD/tmiyamoto/script/for_Kconf_gfiz

mkdir Kconf_${SET}_${CONF}.wdir
cd Kconf_${SET}_${CONF}.wdir
ln -s $KONV96 konv96
ln -s $CoulombGF CoulombGaugeFixingSU3_4D_96x13

sh $SCRIPT/make_Kconfs_bunkatsujob.sh \
    K_${SET}_${CONF}_div \
    Kconf_${SET}_${CONF}_div.out \
    $IDIR $ODIR $CONF $SET > ./Kconf_${SET}_${CONF}_div.job

for TXYZ in t x y z
do
    sh $SCRIPT/make_Kconfs_GFjob.sh \
	Kgf_${SET}_${CONF}${TXYZ} \
	Kconf_${SET}_${CONF}${TXYZ}_gfix.out \
	$ODIR ${TXYZ}${CONF} $PREC $ITER 0 > Kconf_${SET}_${CONF}${TXYZ}_gfix_0.job

    sh $SCRIPT/make_Kconfs_GFjob.sh \
        Kgf_${SET}_${CONF}${TXYZ} \
        Kconf_${SET}_${CONF}${TXYZ}_gfix.out \
        $ODIR ${TXYZ}${CONF} $PREC $ITER 1 > Kconf_${SET}_${CONF}${TXYZ}_gfix_1.job
done

sh $SCRIPT/make_Kconfs_ketsugoujob.sh \
    K_${SET}_${CONF}_comb \
    Kconf_${SET}_${CONF}_comb.out \
    $IDIR $ODIR $CONF $SET > ./Kconf_${SET}_${CONF}_comb.job

ID_div=`echo sleep 10 | qsub Kconf_${SET}_${CONF}_div.job`
ID_gft0=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}t_gfix_0.job`
ID_gfx0=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}x_gfix_0.job`
ID_gfy0=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}y_gfix_0.job`
ID_gfz0=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}z_gfix_0.job`
ID_gft1=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}t_gfix_1.job`
ID_gfx1=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}x_gfix_1.job`
ID_gfy1=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}y_gfix_1.job`
ID_gfz1=`echo sleep 10 | qsub -W depend=afterok:$ID_div Kconf_${SET}_${CONF}z_gfix_1.job`
ID_cmb=`echo sleep 10 | qsub -W depend=afterok:$ID_gft0:$ID_gfx0:$ID_gfy0:$ID_gfz0:$ID_gft1:$ID_gfx1:$ID_gfy1:$ID_gfz1 Kconf_${SET}_${CONF}_comb.job`

cd ../
