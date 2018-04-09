#!/bin/sh

if test $# -ne 6
then
    echo " usage : sh $0 [job name] [out file name] [idir] [odir] \"confs\" [set]"
    exit 0
fi

JOB_NAME=$1
OUT_FILE_NAME=$2
SET_INPUT_PATH=$3
SET_OUTPUT_PATH=$4
SET_CONFS=$5
SET_SET=$6

echo '#PBS -S /bin/bash'
echo "#PBS -N $JOB_NAME"
echo '#PBS -A HALQCD'
echo '#PBS -q comq'
echo '#PBS -l select=1:ncpus=16:mpiprocs=1:ompthreads=16'
echo '#PBS -l walltime=24:00:00'
echo '#PBS -j oe'
echo "#PBS -o $OUT_FILE_NAME"
echo ''
echo '. /opt/Modules/default/init/bash'
echo ''
echo 'module load gnu/4.4.7 cuda/6.5.14'
echo 'module load mvapich2/1.8.1_gnu_cuda-6.0.37'
echo 'module load boost/1.55_gnu-4.4.7'
echo ''
echo 'cd $PBS_O_WORKDIR'
echo ''
echo 'echo $PBS_NODEFILE'
echo 'cat $PBS_NODEFILE'
echo ''
echo '### This job script is made by T.Miyamoto (2015.06.26)'
echo '###    read bunkatsu MDP -> write DDHMC'
echo '### Using, konv96'
echo ''
echo '#########  Parameter Setting #########'
echo '#                                    #'
echo "INPUT_PATH=$SET_INPUT_PATH"
echo "OUTPUT_PATH=$SET_OUTPUT_PATH"
echo ''
echo "CALC_CONF=\"$SET_CONFS\""
echo "SET=$SET_SET"
echo ''
echo 'XYZSIZE=96'
echo 'TSIZE=96'
echo ''
echo '## You need change format "config" ! #'
echo '#                                    #'
echo '######################################'
echo 'DIV_TSIZE=$((TSIZE/8+1))'
echo ''
echo 'for CONF in $CALC_CONF'
echo 'do'
echo ''
echo 'TIME=`date`'
echo 'echo " @@@ MDP->DDHMC START   : conf=$CONF : $TIME"'
echo ''
echo '#####################   bunkatsu MDP -> DDHMC   ######################'
echo '#                                                                    #'
echo 'for DIRECTION in t x y z'
echo 'do'
echo '    mkdir -p $OUTPUT_PATH/Kconfs_gfixed/for_LDDHMC/K8x8x4x4/run_${SET}/${DIRECTION}config_${SET}_${CONF}'
echo '    ./konv96 \'
echo '	-sites $XYZSIZE $XYZSIZE $XYZSIZE $TSIZE \'
echo '	-mdp::bunkatsu_read 8 $OUTPUT_PATH/tmp.mdp/${DIRECTION}${CONF}_t0?.mdp.gfix \'
echo '	-ddhmc::read_profile $INPUT_PATH/$CONF/config.x0y0z0t0 \'
echo '	-ddhmc::NDIM 8 8 4 4 \'
echo '	-default_plaq read_from mdp  \'
echo '	-default_plaq print          \'
echo '	-default_plaq write_to ddhmc \'
echo '	-ddhmc::write $OUTPUT_PATH/Kconfs_gfixed/for_LDDHMC/K8x8x4x4/run_${SET}/${DIRECTION}config_${SET}_${CONF}/${DIRECTION}config_${SET}_${CONF}'
echo '        case "$?" in'
echo '            0) echo "good"; rm ${OUTPUT_PATH}/tmp.mdp/${DIRECTION}${CONF}_t*;;'
echo '            *) echo "@@@@@@ FATAL ERROR !";;'
echo '	esac'
echo 'done'
echo '#                                                                    #'
echo '######################################################################'
echo ''
echo 'TIME=`date`'
echo 'echo " @@@ MDP->DDHMC END     : conf=$CONF : $TIME"'
echo ''
echo 'done'