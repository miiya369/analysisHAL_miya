#!/bin/sh

if [ $# -ne 3 ]; then
    echo "usage: sh `basename $0` [src T] [conf] [txyz]"
    exit 1
fi

SRCT=$1
conf=$2
txyz=$3

srcT=`expr $SRCT + 0`

srcT_F=$((srcT + 1))

CONF=`printf %06d $conf`

echo '// main_quark input file  (config sector)'
echo "$CONF $CONF 10  ! initial traj, final traj, traj skip step"
echo './ConfU.ishikawa.gfix'
echo "${txyz}config_i"
echo '//  hadron/quark save path'
echo "./results.ishikawa.run_i.t_$SRCT.rot_$txyz"
echo '//  (wave function sector)'
echo '3                     ! number of wavefunction'
echo '//  wave function #1 (local)'
echo '1                     ! source smear type: 1=loc, 2=exp, 3=wall, 4=poly'
echo "1 1 1 $srcT_F          ! source center :(x,y,z,t), x=1..NTX, y=1..NTY, z=1..NTZ, t=1..NTT"
echo '//  wave function #2(up/down) (exp)'
echo '3                     ! source smear type: 1=loc, 2=exp, 3=wall, 4=poly'
echo "1 1 1 $srcT_F          ! source center :(x,y,z,t), x=1..NTX, y=1..NTY, z=1..NTZ, t=1..NTT"
echo '//  wave function #3(strange) (exp)'
echo '3                     ! source smear type: 1=loc, 2=exp, 3=wall, 4=poly'
echo "1 1 1 $srcT_F          ! source center :(x,y,z,t), x=1..NTX, y=1..NTY, z=1..NTZ, t=1..NTT"
echo '//  (quark sector)'
echo '2   0.0  1            ! number of quark to be solved [stout_alph  stout_nstep]'
echo '//  quark #1 (up/down)'
echo '0.126117d0 1.11000d0  ! kappa &  csw'
echo '1.0d-06    1000       ! solver tol & maxiter'
echo '1.245d0    1    14    ! SSOR-Omega & SSOR NMR(block) & SAP NSAP'
echo '//  quark #2 (strange)'
echo '0.124790d0 1.11000d0  ! kappa &  csw'
echo '1.0d-12    1000       ! solver tol & maxiter'
echo '1.245d0    1     3    ! SSOR-Omega & SSOR NMR(block) & SAP NSAP'
