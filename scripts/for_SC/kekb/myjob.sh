#!/bin/sh

myname=`whoami`
tmp_dir=`pwd`
bllq > ${tmp_dir}/tmp_queue.txt

s32R=0; s32I=0; s32NQ=0 
a32R=0; a32I=0; a32NQ=0
b32R=0; b32I=0; b32NQ=0
c32R=0; c32I=0; c32NQ=0
d32R=0; d32I=0; d32NQ=0
e32R=0; e32I=0; e32NQ=0
f32R=0; f32I=0; f32NQ=0

while read line
do 
    case "$line" in
	*${myname}*R*32s*) s32R=$((s32R+1));;
	*${myname}*I*32s*) s32I=$((s32I+1));;
	*${myname}*NQ*32s*) s32NQ=$((s32NQ+1));;
        *${myname}*R*32a*) a32R=$((a32R+1));;
	*${myname}*I*32a*) a32I=$((a32I+1));;
	*${myname}*NQ*32a*) a32NQ=$((a32NQ+1));;
        *${myname}*R*32b*) b32R=$((b32R+1));;
	*${myname}*I*32b*) b32I=$((b32I+1));;
	*${myname}*NQ*32b*) b32NQ=$((b32NQ+1));;
        *${myname}*R*32c*) c32R=$((c32R+1));;
	*${myname}*I*32c*) c32I=$((c32I+1));;
	*${myname}*NQ*32c*) c32NQ=$((c32NQ+1));;
        *${myname}*R*32d*) d32R=$((d32R+1));;
	*${myname}*I*32d*) d32I=$((d32I+1));;
	*${myname}*NQ*32d*) d32NQ=$((d32NQ+1));;
        *${myname}*R*32e*) e32R=$((e32R+1));;
	*${myname}*I*32e*) e32I=$((e32I+1));;
	*${myname}*NQ*32e*) e32NQ=$((e32NQ+1));;
        *${myname}*R*32f*) f32R=$((f32R+1));;
	*${myname}*I*32f*) f32I=$((f32I+1));;
	*${myname}*NQ*32f*) f32NQ=$((f32NQ+1));;
    esac
done < ${tmp_dir}/tmp_queue.txt

echo
echo " class| R | I | NQ"
echo " -----+---+---+----"
echo "  32s | $s32R | $s32I | $s32NQ"
echo " -----+---+---+----"
echo "  32a | $a32R | $a32I | $a32NQ"
echo "  32b | $b32R | $b32I | $b32NQ"
echo "  32c | $c32R | $c32I | $c32NQ"
echo " -----+---+---+----"
echo "  32d | $d32R | $d32I | $d32NQ"
echo "  32e | $e32R | $e32I | $e32NQ"
echo "  32f | $f32R | $f32I | $f32NQ"
echo " -----+---+---+----"
echo

rm ${tmp_dir}/tmp_queue.txt
