#!/bin/sh

counter=0

# For S = -1 Octet baryons

#channel=(Pro-Lam Pro-SigZ Neu-SigP)
#prop=(O01O08 O01O04 O02O03)
#
#for snk in `seq 0 1 2`; do
#    for src in `seq 0 1 2`; do
#        snksrc_name="${channel[$snk]}_${channel[$src]}"
#        snksrc_prop="${prop[$snk]}_${prop[$src]}"
#        case $snk in
#            0) had1="proton"; had2="Lambda";;
#            1) had1="proton"; had2="Sigma";;
#            2) had1="proton"; had2="Sigma";;
#        esac
#        num=`printf "%02d" $counter`
#        echo "else if (channel_name == \"$snksrc_name\" ||"
#	echo "         channel_name == \"BBwave.dir.S1.$num\") {"
#        echo "   name      = \"$snksrc_name\";"
#        echo "   directory = \"BBwave.dir.S1.$num\";"
#	echo "   prop_name = \"$snksrc_prop\";"
#	echo "   prop_name += \"_\";"
#	echo "   prop_name += analysis::data_list[SNK_RELA];"
#	echo "   prop_name += \"_\";"
#	echo "   prop_name += analysis::data_list[SRC_RELA];"
#	echo "   OperType  = \"NUC\";"
#        echo "   hadron1.set(\"$had1\");"
#        echo "   hadron2.set(\"$had2\");"
#	echo "   flg_spin_proj = false;"
#        echo "}"
#        counter=$((counter+1))
#    done
#done

# For S = -2 Octet baryons

#channel=(Lam-Lam Pro-XiM Neu-XiZ SigP-SigM SigZ-SigZ SigZ-Lam)
#prop=(O08O08 O01O07 O02O06 O03O05 O04O04 O04O08)
#
#for snk in `seq 0 1 5`; do
#    for src in `seq 0 1 5`; do
#        snksrc_name="${channel[$snk]}_${channel[$src]}"
#        snksrc_prop="${prop[$snk]}_${prop[$src]}"
#        case $snk in
#            0) had1="Lambda"; had2="Lambda";;
#            1) had1="proton"; had2="Xi";;
#            2) had1="proton"; had2="Xi";;
#            3) had1="Sigma"; had2="Sigma";;
#            4) had1="Sigma"; had2="Sigma";;
#	    5) had1="Sigma"; had2="Lambda";;
#        esac
#        num=`printf "%02d" $counter`
#        echo "else if (channel_name == \"$snksrc_name\" ||"
#        echo "         channel_name == \"BBwave.dir.S2.$num\") {"
#        echo "   name      = \"$snksrc_name\";"
#        echo "   directory = \"BBwave.dir.S2.$num\";"
#        echo "   prop_name = \"$snksrc_prop\";"
#        echo "   prop_name += \"_\";"
#        echo "   prop_name += analysis::data_list[SNK_RELA];"
#        echo "   prop_name += \"_\";"
#        echo "   prop_name += analysis::data_list[SRC_RELA];"
#        echo "   OperType  = \"NUC\";"
#        echo "   hadron1.set(\"$had1\");"
#        echo "   hadron2.set(\"$had2\");"
#        echo "   flg_spin_proj = false;"
#        echo "}"
#        counter=$((counter+1))
#    done
#done

# For S = -3 Octet baryons

#channel=(XiZ-Lam XiZ-SigZ XiM-SigP)
#prop=(O06O08 O06O04 O07O03)
#
#for snk in `seq 0 1 2`; do
#    for src in `seq 0 1 2`; do
#        snksrc_name="${channel[$snk]}_${channel[$src]}"
#        snksrc_prop="${prop[$snk]}_${prop[$src]}"
#        case $snk in
#            0) had1="Xi"; had2="Lambda";;
#            1) had1="Xi"; had2="Sigma";;
#            2) had1="Xi"; had2="Sigma";;
#        esac
#        num=`printf "%02d" $counter`
#        echo "else if (channel_name == \"$snksrc_name\" ||"
#        echo "         channel_name == \"BBwave.dir.S3.$num\") {"
#        echo "   name      = \"$snksrc_name\";"
#        echo "   directory = \"BBwave.dir.S3.$num\";"
#        echo "   prop_name = \"$snksrc_prop\";"
#        echo "   prop_name += \"_\";"
#        echo "   prop_name += analysis::data_list[SNK_RELA];"
#        echo "   prop_name += \"_\";"
#        echo "   prop_name += analysis::data_list[SRC_RELA];"
#        echo "   OperType  = \"NUC\";"
#        echo "   hadron1.set(\"$had1\");"
#        echo "   hadron2.set(\"$had2\");"
#        echo "   flg_spin_proj = false;"
#        echo "}"
#        counter=$((counter+1))
#    done
#done

# For S = -1 Isospin projected Octet baryon

channel=(Nuc-Lam12 Nuc-Sig12)

for snk in `seq 0 1 1`; do
    for src in `seq 0 1 1`; do
        snksrc_name="${channel[$snk]}_${channel[$src]}"
        case $snk in
            0) had1="proton"; had2="Lambda";;
            1) had1="proton"; had2="Sigma";;
        esac
        num=`printf "%02d" $counter`
        echo "else if (channel_name == \"$snksrc_name\" ||"
        echo "         channel_name == \"BBwave.prj.S1.$num\") {"
        echo "   name      = \"$snksrc_name\";"
        echo "   directory = \"BBwave.prj.S1.$num\";"
        echo "   OperType  = \"NUC\";"
        echo "   hadron1.set(\"$had1\");"
        echo "   hadron2.set(\"$had2\");"
        echo "   flg_spin_proj = false;"
        echo "}"
        counter=$((counter+1))
    done
done
