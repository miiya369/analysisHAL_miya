#!/bin/sh

if [ $# -ne 3 ]; then
    echo "usage: sh `basename $0` [src T] [conf] [txyz]"
    exit 1
fi

SRCT=$1
conf=$2
txyz=$3

CONF=`printf %06d $conf`

echo '#!/bin/sh -l'
echo ''
echo '#PJM -L rscunit=gwmpc'
echo '#PJM -L rscgrp=batch'
echo '#PJM -L node=256'
echo '#PJM --mpi proc=1024'
echo '#PJM -L elapse=00:30:00'
echo '#PJM --restart'
echo '#PJM -g G16030'
echo "#PJM -o run.out.run_i.t_$SRCT.$CONF.rot_$txyz"
echo '#PJM -j'
echo '#PJM -S'
echo '#PJM --mail-list takaya.miyamoto@yukawa.kyoto-u.ac.jp'
echo '#PJM -m e'
echo ''
echo '## '
echo ''
echo 'module load sparc'
echo 'module load fftw'
echo 'module load boost'
echo ''
echo 'date'
echo 'echo "START"'
echo 'echo ""'
echo ''
echo '########## Edit ##########'
echo 'export PARALLEL=8'
echo 'export OMP_NUM_THREADS=$PARALLEL'
echo ''
echo 'export XOS_MMM_L_ARENA_FREE=2'
echo 'export FLIB_FASTOMP=true'
echo ''
echo 'export rdir_base="./ConfU.ishikawa.gfix"'
echo "export wdir_base=\"./results.ishikawa.run_i.t_$SRCT.rot_$txyz\""
echo ''
echo "export param_input=\"./ParamsHadron.dir/ParamsHadron.t_$SRCT.$CONF.rot_$txyz\""
echo ''
echo '#ln -s ./ConfU $rdir_base'
echo ''
echo '#if [ ! -e $wdir_base/$traj ]; then'
echo '#   mkdir -p $wdir_base/$traj'
echo '#fi'
echo 'if [ ! -e ./results ]; then'
echo '   mkdir -p ./results'
echo 'fi'
echo ''
echo ''
echo 'mpirun \'
echo '    -mca coll_base_reduce_commute_safe 1 \'
echo '    -mca coll_tuned_bcast_same_count 1 \'
echo '    lpgparm -s 4MB -d 4MB -h 4MB -t 4MB -p 4MB \'
echo '    ./main_hadron $param_input'
echo ''
echo 'calc_stat=$?'
echo ''
echo ''
echo "##mv quark_solver_log.${txyz}config_i.$CONF \\"
echo "##   quark_solver_log.run_i.t_$SRCT.$CONF.rot_$txyz"
echo ''
echo ''
echo '##########################'
echo ''
echo 'echo ""'
echo 'echo "END"'
echo 'date'
echo ''
echo 'exit $calc_stat'
