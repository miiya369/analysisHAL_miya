#!/bin/sh

base=/work/HALQCD/tmiyamoto/for_copy/run-h
dire=Kconfs_gfixed/for_LDDHMC/K8x8x4x4/run_h

mrgn="`seq 955 10 985`"
vol1="`seq 995 10 1155`"
vol2="`seq 1165 10 1325`"
vol3="`seq 1335 10 1445`"
vol4="`seq 1455 10 1565`"
vol5="`seq 1575 10 1685`"
vol6="`seq 1695 10 1805`"
vol7="`seq 1815 10 1925`"
vol8="`seq 1935 10 2045`"

Dvol1=Kconfs_gfixed_vol1_68confs_0995-1155
Dvol2=Kconfs_gfixed_vol2_68confs_1165-1325
Dvol3=Kconfs_gfixed_vol3_48confs_1335-1445
Dvol4=Kconfs_gfixed_vol4_48confs_1455-1565
Dvol5=Kconfs_gfixed_vol5_48confs_1575-1685
Dvol6=Kconfs_gfixed_vol6_48confs_1695-1805
Dvol7=Kconfs_gfixed_vol7_48confs_1815-1925
Dvol8=Kconfs_gfixed_vol8_48confs_1935-2045

for conf in $vol1
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol1/$dire/${txyz}config_h_$CONF |wc -w`   
      if test $check -ne 1024
      then
	 echo "ERROR"
      else
	 echo "OK"
      fi
   done
done

for conf in $vol2
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol2/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done

for conf in $vol3
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol3/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done

for conf in $vol4
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol4/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done

for conf in $vol5
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol5/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done

for conf in $vol6
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol6/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done

for conf in $vol7
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol7/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done

for conf in $vol8
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U $base/$Dvol8/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done

for conf in $mrgn
do
   CONF=`printf %06d $conf`
   for txyz in t x y z
   do
      check=`ls -U /work/HALQCD/tmiyamoto/$dire/${txyz}config_h_$CONF |wc -w`
      if test $check -ne 1024
      then
         echo "ERROR"
      else
         echo "OK"
      fi
   done
done