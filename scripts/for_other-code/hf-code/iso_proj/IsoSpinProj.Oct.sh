#!/bin/sh

if [ $# -ne 3 ]; then
    echo "usage: $0 [conf list] [input directory] [output directory]"
    exit 1
fi

######### ! CHANGE HERE ! #########
t_min=12
t_max=15

snk_rela="CG05"
src_rela="CG05"
TXYZshift=A16.000.000.000

# Set YES or NO
Proj_S1_Iso=NO
Proj_S2_Iso=NO
Proj_S3_Iso=NO

AddWave=$HOME/analysis_code-set/main/indep_code/WaveAdd_compress
################################### 

ConfList=$1

ibase=$2
obase=$3

echo "Setup:"
echo " - Current    = `pwd`"
echo " - Conf List  = $ConfList"
echo " - Input dir  = $ibase"
echo " - Output dir = $obase"
echo " - time       = $t_min - $t_max"
echo

if [ $Proj_S1_Iso == "YES" ]; then
    echo
    echo "IsoSpin Projection (S=-1) START: `date`"
    
    sqrt13=`echo "scale=16; sqrt(1/3)" | bc`
    sqrt23=`echo "scale=16; sqrt(2/3)" | bc`
    sqrt19=`echo "scale=16; sqrt(1/9)" | bc`
    sqrt29=`echo "scale=16; sqrt(2/9)" | bc`
    sqrt49=`echo "scale=16; sqrt(4/9)" | bc`
    
    for CONF in `cat $ConfList`
    do
	mkdir -p $obase/BBwave.iso.S1.00/$CONF # N Lambda (I=1/2) - N Lambda (I=1/2)
	mkdir -p $obase/BBwave.iso.S1.01/$CONF # N Lambda (I=1/2) - N  Sigma (I=1/2)
	mkdir -p $obase/BBwave.iso.S1.02/$CONF # N  Sigma (I=1/2) - N Lambda (I=1/2)
	mkdir -p $obase/BBwave.iso.S1.03/$CONF # N  Sigma (I=1/2) - N  Sigma (I=1/2)

	mkdir -p $obase/BBwave.iso.S1.04/$CONF # N  Sigma (I=3/2) - N  Sigma (I=3/2)

	mkdir -p $obase/BBwave.iso.S1.05/$CONF # N Lambda (I=1/2) - N  Sigma (I=3/2) # For Debug
	mkdir -p $obase/BBwave.iso.S1.06/$CONF # N  Sigma (I=1/2) - N  Sigma (I=3/2)
	mkdir -p $obase/BBwave.iso.S1.07/$CONF # N  Sigma (I=3/2) - N  Sigma (I=1/2)
	mkdir -p $obase/BBwave.iso.S1.08/$CONF # N  Sigma (I=3/2) - N Lambda (I=1/2)
	
	for Time in `seq $t_min 1 $t_max`
	do
            TIME=`printf %03d $Time`
	
	    for SRC_RELA in $src_rela
	    do
		for SNK_RELA in $snk_rela
		do
		    BASE=$CONF/NBSwave.+$TIME+$TXYZshift.$CONF.NUC_$SNK_RELA.NUC_$SRC_RELA
		    
		    for i in `seq 0 1 8`; do
			Imat[$i]=$ibase/BBwave.dir.S1.0$i/$BASE
		    done
		    
		    for i in `seq 0 1 8`; do
                        Omat[$i]=$obase/BBwave.iso.S1.0$i/$BASE
		    done
		    
		    $AddWave ${Omat[0]} 1 ${Imat[0]}
		    $AddWave ${Omat[1]} $sqrt13 ${Imat[1]} -$sqrt23 ${Imat[2]}
                    $AddWave ${Omat[2]} $sqrt13 ${Imat[3]} -$sqrt23 ${Imat[6]}
		    $AddWave ${Omat[3]} $sqrt19 ${Imat[4]} -$sqrt29 ${Imat[5]} -$sqrt29 ${Imat[7]} $sqrt49 ${Imat[8]}
		    
		    $AddWave ${Omat[4]} $sqrt19 ${Imat[4]}  $sqrt29 ${Imat[5]}  $sqrt29 ${Imat[7]} $sqrt49 ${Imat[8]}
		    
		    $AddWave ${Omat[5]} $sqrt23 ${Imat[1]}  $sqrt13 ${Imat[2]} # For Debug
                    $AddWave ${Omat[6]} $sqrt29 ${Imat[4]}  $sqrt19 ${Imat[5]} -$sqrt49 ${Imat[7]} -$sqrt29 ${Imat[8]}
                    $AddWave ${Omat[7]} $sqrt29 ${Imat[4]} -$sqrt49 ${Imat[5]}  $sqrt19 ${Imat[7]} -$sqrt29 ${Imat[8]}
		    $AddWave ${Omat[8]} $sqrt23 ${Imat[3]}  $sqrt13 ${Imat[6]}
		done
	    done
	done
	echo "End: $CONF: `date`"
    done
    
    echo "IsoSpin Projection (S=-1)  END : `date`"
fi

if [ $Proj_S2_Iso == "YES" ]; then
    echo
    echo "IsoSpin Projection (S=-2) START: `date`"

    sqrt12=`echo "scale=16; sqrt(1/2)" | bc`
    sqrt13=`echo "scale=16; sqrt(1/3)" | bc`
    sqrt23=`echo "scale=16; sqrt(2/3)" | bc`
    sqrt43=`echo "scale=16; sqrt(4/3)" | bc`
    sqrt16=`echo "scale=16; sqrt(1/6)" | bc`
    sqrt46=`echo "scale=16; sqrt(4/6)" | bc`
    sqrt29=`echo "scale=16; sqrt(2/9)" | bc`
    sqrt89=`echo "scale=16; sqrt(8/9)" | bc`
    over12=`echo "scale=16; (1/2)" | bc`
    over13=`echo "scale=16; (1/3)" | bc`
    over23=`echo "scale=16; (2/3)" | bc`
    over43=`echo "scale=16; (4/3)" | bc`

    for CONF in `cat $ConfList`
    do
        mkdir -p $obase/BBwave.iso.S2.00/$CONF # Lambda Lambda (I=0) - Lambda Lambda (I=0)
        mkdir -p $obase/BBwave.iso.S2.01/$CONF # Lambda Lambda (I=0) -  Sigma  Sigma (I=0)
        mkdir -p $obase/BBwave.iso.S2.02/$CONF # Lambda Lambda (I=0) -      N     Xi (I=0)
        mkdir -p $obase/BBwave.iso.S2.03/$CONF #  Sigma  Sigma (I=0) - Lambda Lambda (I=0)
        mkdir -p $obase/BBwave.iso.S2.04/$CONF #  Sigma  Sigma (I=0) -  Sigma  Sigma (I=0)
        mkdir -p $obase/BBwave.iso.S2.05/$CONF #  Sigma  Sigma (I=0) -      N     Xi (I=0)
        mkdir -p $obase/BBwave.iso.S2.06/$CONF #      N     Xi (I=0) - Lambda Lambda (I=0)
        mkdir -p $obase/BBwave.iso.S2.07/$CONF #      N     Xi (I=0) -  Sigma  Sigma (I=0)
        mkdir -p $obase/BBwave.iso.S2.08/$CONF #      N     Xi (I=0) -      N     Xi (I=0)

        mkdir -p $obase/BBwave.iso.S2.09/$CONF #      N     Xi (I=1) -      N     Xi (I=1)
        mkdir -p $obase/BBwave.iso.S2.10/$CONF #      N     Xi (I=1) -  Sigma Lambda (I=1)
        mkdir -p $obase/BBwave.iso.S2.11/$CONF #  Sigma Lambda (I=1) -      N     Xi (I=1)
        mkdir -p $obase/BBwave.iso.S2.12/$CONF #  Sigma Lambda (I=1) -  Sigma Lambda (I=1)

        mkdir -p $obase/BBwave.iso.S2.13/$CONF #  Sigma  Sigma (I=2) -  Sigma  Sigma (I=2)

        mkdir -p $obase/BBwave.iso.S2.14/$CONF #  Sigma  Sigma (I=0) -  Sigma  Sigma (I=2) # For Debug
        mkdir -p $obase/BBwave.iso.S2.15/$CONF #      N     Xi (I=0) -      N     Xi (I=1)
        mkdir -p $obase/BBwave.iso.S2.16/$CONF # Lambda Lambda (I=0) -  Sigma  Sigma (I=2)
        mkdir -p $obase/BBwave.iso.S2.17/$CONF #  Sigma Lambda (I=1) -      N     Xi (I=0)

        for Time in `seq $t_min 1 $t_max`
        do
            TIME=`printf %03d $Time`

            for SRC_RELA in $src_rela
            do
                for SNK_RELA in $snk_rela
                do
                    BASE=$CONF/NBSwave.+$TIME+$TXYZshift.$CONF.NUC_$SNK_RELA.NUC_$SRC_RELA

                    for i in `seq 0 1 35`; do
                        Imat[$i]=$ibase/BBwave.dir.S2.`printf %02d $i`/$BASE
                    done

                    for i in `seq 0 1 17`; do
                        Omat[$i]=$obase/BBwave.iso.S2.`printf %02d $i`/$BASE
                    done

                    $AddWave ${Omat[0]} 1 ${Imat[0]}
                    $AddWave ${Omat[1]} $sqrt43 ${Imat[ 3]} -$sqrt13 ${Imat[ 4]}
                    $AddWave ${Omat[2]} $sqrt12 ${Imat[ 1]} -$sqrt12 ${Imat[ 2]}
                    $AddWave ${Omat[3]} $sqrt43 ${Imat[18]} -$sqrt13 ${Imat[24]}
                    $AddWave ${Omat[4]} $over43 ${Imat[21]} -$over23 ${Imat[22]} -$over23 ${Imat[27]} $over13 ${Imat[28]}
                    $AddWave ${Omat[5]} $sqrt46 ${Imat[19]} -$sqrt46 ${Imat[20]} -$sqrt16 ${Imat[25]} $sqrt16 ${Imat[26]}
                    $AddWave ${Omat[6]} $sqrt12 ${Imat[ 6]} -$sqrt12 ${Imat[12]}
                    $AddWave ${Omat[7]} $sqrt46 ${Imat[ 9]} -$sqrt16 ${Imat[19]} -$sqrt46 ${Imat[15]} $sqrt16 ${Imat[16]}
                    $AddWave ${Omat[8]} $over12 ${Imat[ 7]} -$over12 ${Imat[ 8]} -$over12 ${Imat[13]} $over12 ${Imat[14]}

                    $AddWave ${Omat[ 9]} $over12 ${Imat[ 7]} $over12 ${Imat[ 8]} $over12 ${Imat[13]} $over12 ${Imat[14]}
                    $AddWave ${Omat[10]} $sqrt12 ${Imat[11]} $sqrt12 ${Imat[17]}
                    $AddWave ${Omat[11]} $sqrt12 ${Imat[31]} $sqrt12 ${Imat[32]}
                    $AddWave ${Omat[12]} 1 ${Imat[35]}
		    
                    $AddWave ${Omat[13]} $over23 ${Imat[21]} $over23 ${Imat[22]} $over23 ${Imat[27]} $over23 ${Imat[28]}
		    
		    # For Debug
                    $AddWave ${Omat[14]} $sqrt89 ${Imat[21]}  $sqrt89 ${Imat[22]} -$sqrt29 ${Imat[27]} -$sqrt29 ${Imat[28]}
                    $AddWave ${Omat[15]} $over12 ${Imat[ 7]}  $over12 ${Imat[ 8]} -$over12 ${Imat[13]} -$over12 ${Imat[14]}
                    $AddWave ${Omat[16]} $sqrt23 ${Imat[ 3]}  $sqrt23 ${Imat[ 4]}
                    $AddWave ${Omat[17]} $sqrt12 ${Imat[31]} -$sqrt12 ${Imat[32]}
                done
            done
        done
        echo "End: $CONF: `date`"
    done

    echo "IsoSpin Projection (S=-2)  END : `date`"
fi

if [ $Proj_S3_Iso == "YES" ]; then
    echo
    echo "IsoSpin Projection (S=-3) START: `date`"

    sqrt13=`echo "scale=16; sqrt(1/3)" | bc`
    sqrt23=`echo "scale=16; sqrt(2/3)" | bc`
    sqrt19=`echo "scale=16; sqrt(1/9)" | bc`
    sqrt29=`echo "scale=16; sqrt(2/9)" | bc`
    sqrt49=`echo "scale=16; sqrt(4/9)" | bc`

    for CONF in `cat $ConfList`
    do
        mkdir -p $obase/BBwave.iso.S3.00/$CONF # Xi Lambda (I=1/2) - Xi Lambda (I=1/2)
        mkdir -p $obase/BBwave.iso.S3.01/$CONF # Xi Lambda (I=1/2) - Xi  Sigma (I=1/2)
        mkdir -p $obase/BBwave.iso.S3.02/$CONF # Xi  Sigma (I=1/2) - Xi Lambda (I=1/2)
        mkdir -p $obase/BBwave.iso.S3.03/$CONF # Xi  Sigma (I=1/2) - Xi  Sigma (I=1/2)

        mkdir -p $obase/BBwave.iso.S3.04/$CONF # Xi  Sigma (I=3/2) - Xi  Sigma (I=3/2)

        mkdir -p $obase/BBwave.iso.S3.05/$CONF # Xi Lambda (I=1/2) - Xi  Sigma (I=3/2) # For Debug
        mkdir -p $obase/BBwave.iso.S3.06/$CONF # Xi  Sigma (I=1/2) - Xi  Sigma (I=3/2)
        mkdir -p $obase/BBwave.iso.S3.07/$CONF # Xi  Sigma (I=3/2) - Xi  Sigma (I=1/2)
        mkdir -p $obase/BBwave.iso.S3.08/$CONF # Xi  Sigma (I=3/2) - Xi Lambda (I=1/2)

        for Time in `seq $t_min 1 $t_max`
        do
            TIME=`printf %03d $Time`

            for SRC_RELA in $src_rela
            do
                for SNK_RELA in $snk_rela
                do
                    BASE=$CONF/NBSwave.+$TIME+$TXYZshift.$CONF.NUC_$SNK_RELA.NUC_$SRC_RELA

                    for i in `seq 0 1 8`; do
                        Imat[$i]=$ibase/BBwave.dir.S3.0$i/$BASE
                    done

                    for i in `seq 0 1 8`; do
                        Omat[$i]=$obase/BBwave.iso.S3.0$i/$BASE
                    done

                    $AddWave ${Omat[0]} 1 ${Imat[0]}
                    $AddWave ${Omat[1]} $sqrt13 ${Imat[1]} -$sqrt23 ${Imat[2]}
                    $AddWave ${Omat[2]} $sqrt13 ${Imat[3]} -$sqrt23 ${Imat[6]}
                    $AddWave ${Omat[3]} $sqrt19 ${Imat[4]} -$sqrt29 ${Imat[5]} -$sqrt29 ${Imat[7]} $sqrt49 ${Imat[8]}
                    
		    $AddWave ${Omat[4]} $sqrt19 ${Imat[4]}  $sqrt29 ${Imat[5]}  $sqrt29 ${Imat[7]} $sqrt49 ${Imat[8]}

                    $AddWave ${Omat[5]} $sqrt23 ${Imat[1]}  $sqrt13 ${Imat[2]} # For Debug
                    $AddWave ${Omat[6]} $sqrt29 ${Imat[4]}  $sqrt19 ${Imat[5]} -$sqrt49 ${Imat[7]} -$sqrt29 ${Imat[8]}
                    $AddWave ${Omat[7]} $sqrt29 ${Imat[4]} -$sqrt49 ${Imat[5]}  $sqrt19 ${Imat[7]} -$sqrt29 ${Imat[8]}
                    $AddWave ${Omat[8]} $sqrt23 ${Imat[3]}  $sqrt13 ${Imat[6]}
                done
            done
        done
        echo "End: $CONF: `date`"
    done

    echo "IsoSpin Projection (S=-3)  END : `date`"
fi

################################################################
########################### For Test ###########################

##Pmat=("1"       "0"        "0" \
##      "0" "Sq(1/3)" "-Sq(2/3)" \
##      "0" "Sq(2/3)"  "Sq(1/3)")
#
#Pmat=("1"       "0"        "0"        "0"       "0"        "0" "0" \
#      "0" "Sq(1/2)" "-Sq(1/2)"        "0"       "0"        "0" "0" \
#      "0" "Sq(1/2)"  "Sq(1/2)"        "0"       "0"        "0" "0" \
#      "0"       "0"        "0" "Sq(1/3)" "-Sq(1/3)"  "Sq(1/3)" "0" \
#      "0"       "0"        "0" "Sq(1/2)"        "0" "-Sq(1/2)" "0" \
#      "0"       "0"        "0" "Sq(1/6)"  "Sq(2/3)"  "Sq(1/6)" "0" \
#      "0"       "0"        "0"        "0"       "0"        "0" "1")
#
##N=3
#N=7
#
#NN=$((N*N))
#for i in `seq 0 1 $((NN-1))`; do
#    Imat[$i]="d$i"
#done
#
#for i in `seq 0 1 $((N-1))`; do
#    for j in `seq 0 1 $((N-1))`; do
#	for k in `seq 0 1 $((N-1))`; do
#	    for l in `seq 0 1 $((N-1))`; do
#		if [ ${Pmat[$((k+N*i))]} != "0" -a ${Pmat[$((l+N*j))]} != "0" ]; then 
#		    Omat[$((j+N*i))]="${Omat[$((j+N*i))]} + ${Pmat[$((k+N*i))]}x${Imat[$((l+N*k))]}x${Pmat[$((l+N*j))]}"
#		fi
#	    done
#	done
#    done
#done
#
#for i in `seq 0 1 $((NN-1))`; do
#    echo ${Omat[$i]}
#done
