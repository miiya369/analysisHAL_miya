#!/bin/sh

if [ $# -ne 1 ]; then
    echo "usage: sh `basename $0` [Tsrc % 24 = ?]"
    exit 1
fi

TsrcMod24=$1
TSRCMOD24=`printf %03d $TsrcMod24`

Base1=/data/G15023/gongyo/confs/kei/K8x8x4x4/run-i/traj_mod10_05/gfix_cond_1.0e-14
Base2=../../Jobs.dir.run.v04.for_FX100.run-i.doi.2016-04-17/Jobs.t_mod24_$TSRCMOD24
Base3=../../ParamsHadron.dir.ud_e-06.s_e-12.v03.for_FX100.run-i.doi.2015-06-14/Params.t_mod24_$TSRCMOD24
Base4=../../../hadron_force/hadron_force.2016-07-06/tests/test.loop.lst/loop.lst
Base5=../../../ishikawa_solver/LDDHMC_ver1.3_K0.52_ovlpcomm1.20_stout-r309.hal.v13/Measurement/HadronBlock/main_hadron
Base6=../../../hadron_force/hadron_force.2016-07-06/tests/test.hadron_force.2b/work

for Tsrc in `seq $TsrcMod24 24 95`; do
    TSRC=`printf %03d $Tsrc`
    
    mkdir run.t_$TSRC.rot_tot
    cd run.t_$TSRC.rot_tot
    ln -s $Base1 ConfU
    ln -s ConfU ConfU.ishikawa.gfix
    ln -s $Base2 Jobs.dir
    ln -s $Base3 ParamsHadron.dir
    ln -s $Base4 loop.lst
    ln -s $Base5 main_hadron
    ln -s $Base6 work
    ln -s results results.t_$TSRC
    mkdir logs.t_$TSRC
    cd ../
done
