#!/bin/sh

OK=0
NO=0

echo " Expect : #.file = 1024"

count=`ls|wc -w`
conf=`ls`

echo -n " Start"
for i in `seq 1 1 $count`
do
   if test $((i % 5)) -eq 0
   then echo -n "+"
   else echo -n "-"
   fi
done
echo "End"
echo -n "      "

for CONF in $conf
do 
   tmp=`ls -U $CONF | wc -w`
   if test $tmp -eq 1024
   then
      OK=$((OK+1))
   else
      NO=$((NO+1))
   fi
   echo -n "*"
done

echo ""
echo " OK = `printf %03d $OK` / `printf %03d $count`"
echo " NO = `printf %03d $NO` / `printf %03d $count`"