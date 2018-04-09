#!/bin/bash

if [ $# -ne 3 ]; then
    echo "usage: $0 [conf list] [input directory] [output directory]"
    exit 1
fi

######### ! CHANGE HERE ! #########
t_min=6
t_max=12

snk_rela="NR"
src_rela="NR"
TXYZshift=A00.000.000.000

AddWave=$HOME/analysisHAL_miya/main/indep_code/WaveAdd_not-compress

sub_dir="spin1_ave"
################################### 

ConfList=$1

ibase=$2
obase=$3

echo "Setup:"
echo " - Current    = `pwd`"
echo " - Conf List  = $ConfList"
echo " - Input dir  = $ibase"
echo " - Output dir = $obase"
echo " - Time       = $t_min - $t_max"
echo " - Sub dir    = $sub_dir"
echo 

sqrt13=`echo "scale=16; sqrt(1/3)" | bc`
sqrt23=`echo "scale=16; sqrt(2/3)" | bc`
sqrt19=`echo "scale=16; sqrt(1/9)" | bc`
sqrt29=`echo "scale=16; sqrt(2/9)" | bc`
sqrt49=`echo "scale=16; sqrt(4/9)" | bc`

echo "JOB START: `date`"
echo

for CONF in `cat $ConfList`
do
    for SubDir in $sub_dir
    do
	NL___NL__=$obase/IsoProj.NBS_S1.Nuc__Lam__Nuc__Lam_.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NL___NL__
	NL___NS1_=$obase/IsoProj.NBS_S1.Nuc__Lam__Nuc__Sig_.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NL___NS1_
	NL___NSS1=$obase/IsoProj.NBS_S1.Nuc__Lam__Nuc__SigS.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NL___NSS1
	
	NS1__NL__=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__Lam_.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NS1__NL__
	NS1__NS1_=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__Sig_.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NS1__NS1_
	NS1__NSS1=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__SigS.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NS1__NSS1
	
	NSS1_NL__=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__Lam_.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NSS1_NL__
	NSS1_NS1_=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__Sig_.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NSS1_NS1_
	NSS1_NSS1=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__SigS.dir/iso.1p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NSS1_NSS1
	
	NS3__NS3_=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__Sig_.dir/iso.3p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NS3__NS3_
	NS3__NSS3=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__SigS.dir/iso.3p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NS3__NSS3
	
	NSS3_NS3_=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__Sig_.dir/iso.3p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NSS3_NS3_
	NSS3_NSS3=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__SigS.dir/iso.3p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NSS3_NSS3
	
	### For Debug ###
	NL___NS3_=$obase/IsoProj.NBS_S1.Nuc__Lam__Nuc__Sig_.dir/iso.1p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NL___NS3_
	NL___NSS3=$obase/IsoProj.NBS_S1.Nuc__Lam__Nuc__SigS.dir/iso.1p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NL___NSS3
	
	NS3__NL__=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__Lam_.dir/iso.3p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NS3__NL__
	NS3__NS1_=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__Sig_.dir/iso.3p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NS3__NS1_
	NS3__NSS1=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__SigS.dir/iso.3p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NS3__NSS1
	
	NSS3_NL__=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__Lam_.dir/iso.3p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NSS3_NL__
	NSS3_NS1_=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__Sig_.dir/iso.3p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NSS3_NS1_
	NSS3_NSS1=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__SigS.dir/iso.3p2z+1p2.1p2z+1p2/$SubDir/$CONF; mkdir -p $NSS3_NSS1
	
	NS1__NS3_=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__Sig_.dir/iso.1p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NS1__NS3_
	NS1__NSS3=$obase/IsoProj.NBS_S1.Nuc__Sig__Nuc__SigS.dir/iso.1p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NS1__NSS3
	
	NSS1_NS3_=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__Sig_.dir/iso.1p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NSS1_NS3_
	NSS1_NSS3=$obase/IsoProj.NBS_S1.Nuc__SigS_Nuc__SigS.dir/iso.1p2z+1p2.3p2z+1p2/$SubDir/$CONF; mkdir -p $NSS1_NSS3
	#################
	
	for Time in `seq $t_min 1 $t_max`
	do
            TIME=`printf %03d $Time`
	    
	    for SRC_RELA in $src_rela
	    do
		for SNK_RELA in $snk_rela
		do
		    iBASE=NBSwave.+$TIME+$TXYZshift.$CONF.??_$SNK_RELA.??_$SRC_RELA
		    
		    oBASE1=NBSwave.+$TIME+$TXYZshift.$CONF.snk_$SNK_RELA.src_$SRC_RELA
		    oBASE2=NBSwave.+$TIME+$TXYZshift.$CONF.snk_$SNK_RELA.src_$SRC_RELA
		    oBASE3=NBSwave.+$TIME+$TXYZshift.$CONF.snk_$SNK_RELA.src_$SRC_RELA
		    oBASE4=NBSwave.+$TIME+$TXYZshift.$CONF.snk_$SNK_RELA.src_$SRC_RELA
		    
		    pL___pL__=$ibase/Proj.NBS_2Boct.Prot__Lamb__Prot__Lamb_.dir/$SubDir/$CONF/$iBASE
		    pL___pSZ_=$ibase/Proj.NBS_2Boct.Prot__Lamb__Prot__SigZ_.dir/$SubDir/$CONF/$iBASE
		    pL___nSP_=$ibase/Proj.NBS_2Boct.Prot__Lamb__Neut__SigP_.dir/$SubDir/$CONF/$iBASE
		    pL___pSSZ=$ibase/Proj.NBS_ooxdo.Prot__Lamb__SigSZ_Prot_.dir/$SubDir/$CONF/$iBASE
		    pL___nSSP=$ibase/Proj.NBS_ooxdo.Prot__Lamb__SigSP_Neut_.dir/$SubDir/$CONF/$iBASE
		    pSZ__pL__=$ibase/Proj.NBS_2Boct.Prot__SigZ__Prot__Lamb_.dir/$SubDir/$CONF/$iBASE
		    pSZ__pSZ_=$ibase/Proj.NBS_2Boct.Prot__SigZ__Prot__SigZ_.dir/$SubDir/$CONF/$iBASE
		    pSZ__nSP_=$ibase/Proj.NBS_2Boct.Prot__SigZ__Neut__SigP_.dir/$SubDir/$CONF/$iBASE
		    pSZ__pSSZ=$ibase/Proj.NBS_ooxdo.Prot__SigZ__SigSZ_Prot_.dir/$SubDir/$CONF/$iBASE
		    pSZ__nSSP=$ibase/Proj.NBS_ooxdo.Prot__SigZ__SigSP_Neut_.dir/$SubDir/$CONF/$iBASE
		    nSP__pL__=$ibase/Proj.NBS_2Boct.Neut__SigP__Prot__Lamb_.dir/$SubDir/$CONF/$iBASE
		    nSP__pSZ_=$ibase/Proj.NBS_2Boct.Neut__SigP__Prot__SigZ_.dir/$SubDir/$CONF/$iBASE
		    nSP__nSP_=$ibase/Proj.NBS_2Boct.Neut__SigP__Neut__SigP_.dir/$SubDir/$CONF/$iBASE
		    nSP__pSSZ=$ibase/Proj.NBS_ooxdo.Neut__SigP__SigSZ_Prot_.dir/$SubDir/$CONF/$iBASE
		    nSP__nSSP=$ibase/Proj.NBS_ooxdo.Neut__SigP__SigSP_Neut_.dir/$SubDir/$CONF/$iBASE
		    pSSZ_pL__=$ibase/Proj.NBS_doxoo.SigSZ_Prot__Prot__Lamb_.dir/$SubDir/$CONF/$iBASE
		    pSSZ_pSZ_=$ibase/Proj.NBS_doxoo.SigSZ_Prot__Prot__SigZ_.dir/$SubDir/$CONF/$iBASE
		    pSSZ_nSP_=$ibase/Proj.NBS_doxoo.SigSZ_Prot__Neut__SigP_.dir/$SubDir/$CONF/$iBASE
		    pSSZ_pSSZ=$ibase/Proj.NBS_2Bdow.SigSZ_Prot__SigSZ_Prot_.dir/$SubDir/$CONF/$iBASE
		    pSSZ_nSSP=$ibase/Proj.NBS_2Bdow.SigSZ_Prot__SigSP_Neut_.dir/$SubDir/$CONF/$iBASE
		    nSSP_pL__=$ibase/Proj.NBS_doxoo.SigSP_Neut__Prot__Lamb_.dir/$SubDir/$CONF/$iBASE
		    nSSP_pSZ_=$ibase/Proj.NBS_doxoo.SigSP_Neut__Prot__SigZ_.dir/$SubDir/$CONF/$iBASE
		    nSSP_nSP_=$ibase/Proj.NBS_doxoo.SigSP_Neut__Neut__SigP_.dir/$SubDir/$CONF/$iBASE
		    nSSP_pSSZ=$ibase/Proj.NBS_2Bdow.SigSP_Neut__SigSZ_Prot_.dir/$SubDir/$CONF/$iBASE
		    nSSP_nSSP=$ibase/Proj.NBS_2Bdow.SigSP_Neut__SigSP_Neut_.dir/$SubDir/$CONF/$iBASE
		    
		    $AddWave $NL___NL__/$oBASE1      1.0 $pL___pL__
		    $AddWave $NL___NS1_/$oBASE1 -$sqrt23 $pL___nSP_  $sqrt13 $pL___pSZ_
		    $AddWave $NL___NSS1/$oBASE2 -$sqrt23 $pL___nSSP  $sqrt13 $pL___pSSZ
		    
		    $AddWave $NS1__NL__/$oBASE1 -$sqrt23 $nSP__pL__  $sqrt13 $pSZ__pL__
		    $AddWave $NS1__NS1_/$oBASE1  $sqrt49 $nSP__nSP_ -$sqrt29 $nSP__pSZ_ -$sqrt29 $pSZ__nSP_  $sqrt19 $pSZ__pSZ_
		    $AddWave $NS1__NSS1/$oBASE2  $sqrt49 $nSP__nSSP -$sqrt29 $nSP__pSSZ -$sqrt29 $pSZ__nSSP  $sqrt19 $pSZ__pSSZ
		    
		    $AddWave $NSS1_NL__/$oBASE3 -$sqrt23 $nSSP_pL__  $sqrt13 $pSSZ_pL__
		    $AddWave $NSS1_NS1_/$oBASE3  $sqrt49 $nSSP_nSP_ -$sqrt29 $nSSP_pSZ_ -$sqrt29 $pSSZ_nSP_  $sqrt19 $pSSZ_pSZ_
		    $AddWave $NSS1_NSS1/$oBASE4  $sqrt49 $nSSP_nSSP -$sqrt29 $nSSP_pSSZ -$sqrt29 $pSSZ_nSSP  $sqrt19 $pSSZ_pSSZ
		    
		    $AddWave $NS3__NS3_/$oBASE1  $sqrt19 $nSP__nSP_  $sqrt29 $nSP__pSZ_  $sqrt29 $pSZ__nSP_  $sqrt49 $pSZ__pSZ_
		    $AddWave $NS3__NSS3/$oBASE2  $sqrt19 $nSP__nSSP  $sqrt29 $nSP__pSSZ  $sqrt29 $pSZ__nSSP  $sqrt49 $pSZ__pSSZ
		    
		    $AddWave $NSS3_NS3_/$oBASE3  $sqrt19 $nSSP_nSP_  $sqrt29 $nSSP_pSZ_  $sqrt29 $pSSZ_nSP_  $sqrt49 $pSSZ_pSZ_
		    $AddWave $NSS3_NSS3/$oBASE4  $sqrt19 $nSSP_nSSP  $sqrt29 $nSSP_pSSZ  $sqrt29 $pSSZ_nSSP  $sqrt49 $pSSZ_pSSZ
		    
		    ### For Debug ###
		    $AddWave $NL___NS3_/$oBASE1  $sqrt13 $pL___nSP_  $sqrt23 $pL___pSZ_
		    $AddWave $NL___NSS3/$oBASE2  $sqrt13 $pL___nSSP  $sqrt23 $pL___pSSZ
		    
		    $AddWave $NS3__NL__/$oBASE1  $sqrt13 $nSP__pL__  $sqrt23 $pSZ__pL__
		    $AddWave $NS3__NS1_/$oBASE1 -$sqrt29 $nSP__nSP_  $sqrt19 $nSP__pSZ_ -$sqrt49 $pSZ__nSP_  $sqrt29 $pSZ__pSZ_
		    $AddWave $NS3__NSS1/$oBASE2 -$sqrt29 $nSP__nSSP  $sqrt19 $nSP__pSSZ -$sqrt49 $pSZ__nSSP  $sqrt29 $pSZ__pSSZ
		    
		    $AddWave $NSS3_NL__/$oBASE3  $sqrt13 $nSSP_pL__  $sqrt23 $pSSZ_pL__
		    $AddWave $NSS3_NS1_/$oBASE3 -$sqrt29 $nSSP_nSP_  $sqrt19 $nSSP_pSZ_ -$sqrt49 $pSSZ_nSP_  $sqrt29 $pSSZ_pSZ_
		    $AddWave $NSS3_NSS1/$oBASE4 -$sqrt29 $nSSP_nSSP  $sqrt19 $nSSP_pSSZ -$sqrt49 $pSSZ_nSSP  $sqrt29 $pSSZ_pSSZ
		    
		    $AddWave $NS1__NS3_/$oBASE1 -$sqrt29 $nSP__nSP_ -$sqrt49 $nSP__pSZ_  $sqrt19 $pSZ__nSP_  $sqrt29 $pSZ__pSZ_
		    $AddWave $NS1__NSS3/$oBASE2 -$sqrt29 $nSP__nSSP -$sqrt49 $nSP__pSSZ  $sqrt19 $pSZ__nSSP  $sqrt29 $pSZ__pSSZ
		    
		    $AddWave $NSS1_NS3_/$oBASE3 -$sqrt29 $nSSP_nSP_ -$sqrt49 $nSSP_pSZ_  $sqrt19 $pSSZ_nSP_  $sqrt29 $pSSZ_pSZ_
		    $AddWave $NSS1_NSS3/$oBASE4 -$sqrt29 $nSSP_nSSP -$sqrt49 $nSSP_pSSZ  $sqrt19 $pSSZ_nSSP  $sqrt29 $pSSZ_pSSZ
		    #################
		done
	    done
	done
    done
    echo "Conf = $CONF DONE: `date`"
done

echo
echo "JOB  END : `date`"
