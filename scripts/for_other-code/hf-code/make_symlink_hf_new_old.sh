#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: sh `basename $0` [ABSOLUTE path for results of old HF-code)] [Path for symlink dir]"; exit -1
fi

ibase=$1
obase=$2

nbs_dir_old=(BBwave.dir.S0.00 BBwave.dir.S1.00 BBwave.dir.S1.01 BBwave.dir.S1.02 BBwave.dir.S1.03 BBwave.dir.S1.04 BBwave.dir.S1.05 BBwave.dir.S1.06 BBwave.dir.S1.07 BBwave.dir.S1.08 BBwave.dir.S2.00 BBwave.dir.S2.01 BBwave.dir.S2.02 BBwave.dir.S2.03 BBwave.dir.S2.04 BBwave.dir.S2.06 BBwave.dir.S2.07 BBwave.dir.S2.08 BBwave.dir.S2.09 BBwave.dir.S2.10 BBwave.dir.S2.11 BBwave.dir.S2.12 BBwave.dir.S2.13 BBwave.dir.S2.14 BBwave.dir.S2.15 BBwave.dir.S2.16 BBwave.dir.S2.17 BBwave.dir.S2.18 BBwave.dir.S2.19 BBwave.dir.S2.20 BBwave.dir.S2.21 BBwave.dir.S2.22 BBwave.dir.S2.23 BBwave.dir.S2.24 BBwave.dir.S2.25 BBwave.dir.S2.26 BBwave.dir.S2.27 BBwave.dir.S2.28 BBwave.dir.S2.31 BBwave.dir.S2.32 BBwave.dir.S2.33 BBwave.dir.S2.35 BBwave.dir.S3.00 BBwave.dir.S3.01 BBwave.dir.S3.02 BBwave.dir.S3.03 BBwave.dir.S3.04 BBwave.dir.S3.05 BBwave.dir.S3.06 BBwave.dir.S3.07 BBwave.dir.S3.08 BBwave.dir.S4.00 PPwave.dir.S0.00 PPwave.dir.S1.00 PPwave.dir.S2.00 Proj.DDwave.dir.S6.00 Proj.DOwave.dir.S3.35 Proj.DOwave.dir.S5.00)

nbs_dir_new=(NBS_2Boct.Prot__Neut__Prot__Neut_.dir NBS_2Boct.Prot__Lamb__Prot__Lamb_.dir NBS_2Boct.Prot__Lamb__Prot__SigZ_.dir NBS_2Boct.Prot__Lamb__Neut__SigP_.dir NBS_2Boct.Prot__SigZ__Prot__Lamb_.dir NBS_2Boct.Prot__SigZ__Prot__SigZ_.dir NBS_2Boct.Prot__SigZ__Neut__SigP_.dir NBS_2Boct.Neut__SigP__Prot__Lamb_.dir NBS_2Boct.Neut__SigP__Prot__SigZ_.dir NBS_2Boct.Neut__SigP__Neut__SigP_.dir NBS_2Boct.Lamb__Lamb__Lamb__Lamb_.dir NBS_2Boct.Lamb__Lamb__Prot__XiM__.dir NBS_2Boct.Lamb__Lamb__Neut__XiZ__.dir NBS_2Boct.Lamb__Lamb__SigP__SigM_.dir NBS_2Boct.Lamb__Lamb__SigZ__SigZ_.dir NBS_2Boct.Prot__XiM___Lamb__Lamb_.dir NBS_2Boct.Prot__XiM___Prot__XiM__.dir NBS_2Boct.Prot__XiM___Neut__XiZ__.dir NBS_2Boct.Prot__XiM___SigP__SigM_.dir NBS_2Boct.Prot__XiM___SigZ__SigZ_.dir NBS_2Boct.Prot__XiM___SigZ__Lamb_.dir NBS_2Boct.Neut__XiZ___Lamb__Lamb_.dir NBS_2Boct.Neut__XiZ___Prot__XiM__.dir NBS_2Boct.Neut__XiZ___Neut__XiZ__.dir NBS_2Boct.Neut__XiZ___SigP__SigM_.dir NBS_2Boct.Neut__XiZ___SigZ__SigZ_.dir NBS_2Boct.Neut__XiZ___SigZ__Lamb_.dir NBS_2Boct.SigP__SigM__Lamb__Lamb_.dir NBS_2Boct.SigP__SigM__Prot__XiM__.dir NBS_2Boct.SigP__SigM__Neut__XiZ__.dir NBS_2Boct.SigP__SigM__SigP__SigM_.dir NBS_2Boct.SigP__SigM__SigZ__SigZ_.dir NBS_2Boct.SigP__SigM__SigZ__Lamb_.dir NBS_2Boct.SigZ__SigZ__Lamb__Lamb_.dir NBS_2Boct.SigZ__SigZ__Prot__XiM__.dir NBS_2Boct.SigZ__SigZ__Neut__XiZ__.dir NBS_2Boct.SigZ__SigZ__SigP__SigM_.dir NBS_2Boct.SigZ__SigZ__SigZ__SigZ_.dir NBS_2Boct.SigZ__Lamb__Prot__XiM__.dir NBS_2Boct.SigZ__Lamb__Neut__XiZ__.dir NBS_2Boct.SigZ__Lamb__SigP__SigM_.dir NBS_2Boct.SigZ__Lamb__SigZ__Lamb_.dir NBS_2Boct.XiZ___Lamb__XiZ___Lamb_.dir NBS_2Boct.XiZ___Lamb__XiZ___SigZ_.dir NBS_2Boct.XiZ___Lamb__XiM___SigP_.dir NBS_2Boct.XiZ___SigZ__XiZ___Lamb_.dir NBS_2Boct.XiZ___SigZ__XiZ___SigZ_.dir NBS_2Boct.XiZ___SigZ__XiM___SigP_.dir NBS_2Boct.XiM___SigP__XiZ___Lamb_.dir NBS_2Boct.XiM___SigP__XiZ___SigZ_.dir NBS_2Boct.XiM___SigP__XiM___SigP_.dir NBS_2Boct.XiZ___XiM___XiZ___XiM__.dir NBS_2Mpsw.Pion__Pion__Pion__Pion_.dir NBS_2Mpsw.Pion__Kaon__Pion__Kaon_.dir NBS_2Mpsw.Kaon__Kaon__Kaon__Kaon_.dir Proj.NBS_2Bdec.OmgM__OmgM__OmgM__OmgM_.dir Proj.NBS_2Bdow.OmgM__Prot__OmgM__Prot_.dir Proj.NBS_2Bdow.OmgM__XiZ___OmgM__XiZ__.dir)

echo "${#nbs_dir_old[@]} ${#nbs_dir_new[@]}"
if [ ! -d $obase ]; then
    mkdir -p $obase
fi
for i in `seq 0 1 $((${#nbs_dir_old[@]}-1))`; do
    echo "Make symlink for ${nbs_dir_old[$i]} -> ${nbs_dir_new[$i]}"
    ln -s $ibase/${nbs_dir_old[$i]} $obase/${nbs_dir_new[$i]}
done
