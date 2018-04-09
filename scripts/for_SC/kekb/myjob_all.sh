#!/bin/sh

myname=`whoami`
tmp_dir=`pwd`
bllq > ${tmp_dir}/tmp_queue.txt

s32R=0; s32I=0; s32NQ=0; a32R=0; a32I=0; a32NQ=0; b32R=0; b32I=0; b32NQ=0
c32R=0; c32I=0; c32NQ=0; d32R=0; d32I=0; d32NQ=0; e32R=0; e32I=0; e32NQ=0
f32R=0; f32I=0; f32NQ=0; m32R=0; m32I=0; m32NQ=0

a128R=0; a128I=0; a128NQ=0; b128R=0; b128I=0; b128NQ=0; c128R=0; c128I=0; c128NQ=0
d128R=0; d128I=0; d128NQ=0; e128R=0; e128I=0; e128NQ=0; f128R=0; f128I=0; f128NQ=0
g128R=0; g128I=0; g128NQ=0; h128R=0; h128I=0; h128NQ=0; i128R=0; i128I=0; i128NQ=0
j128R=0; j128I=0; j128NQ=0

a512R=0; a512I=0; a512NQ=0; b512R=0; b512I=0; b512NQ=0; c512R=0; c512I=0; c512NQ=0
d512R=0; d512I=0; d512NQ=0; e512R=0; e512I=0; e512NQ=0

a2048R=0; a2048I=0; a2048NQ=0

while read line
do 
    case "$line" in
	*${myname}*R*32s*) s32R=$((s32R+1));;
	*${myname}*I*32s*) s32I=$((s32I+1));;
	*${myname}*NQ*32s*) s32NQ=$((s32NQ+1));;
	*${myname}*R*32m*) m32R=$((m32R+1));;
        *${myname}*I*32m*) m32I=$((m32I+1));;
        *${myname}*NQ*32m*) m32NQ=$((m32NQ+1));;
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
        *${myname}*R*128a*) a128R=$((a128R+1));;
        *${myname}*I*128a*) a128I=$((a128I+1));;
        *${myname}*NQ*128a*) a128NQ=$((a128NQ+1));;
        *${myname}*R*128b*) b128R=$((b128R+1));;
        *${myname}*I*128b*) b128I=$((b128I+1));;
        *${myname}*NQ*128b*) b128NQ=$((b128NQ+1));;
        *${myname}*R*128c*) c128R=$((c128R+1));;
        *${myname}*I*128c*) c128I=$((c128I+1));;
        *${myname}*NQ*128c*) c128NQ=$((c128NQ+1));;
        *${myname}*R*128d*) d128R=$((d128R+1));;
        *${myname}*I*128d*) d128I=$((d128I+1));;
        *${myname}*NQ*128d*) d128NQ=$((d128NQ+1));;
        *${myname}*R*128e*) e128R=$((e128R+1));;
        *${myname}*I*128e*) e128I=$((e128I+1));;
        *${myname}*NQ*128e*) e128NQ=$((e128NQ+1));;
        *${myname}*R*128f*) f128R=$((f128R+1));;
	*${myname}*I*128f*) f128I=$((f128I+1));;
	*${myname}*NQ*128f*) f128NQ=$((f128NQ+1));;
        *${myname}*R*128g*) g128R=$((g128R+1));;
	*${myname}*I*128g*) g128I=$((g128I+1));;
        *${myname}*NQ*128g*) g128NQ=$((g128NQ+1));;
        *${myname}*R*128h*) h128R=$((h128R+1));;
        *${myname}*I*128h*) h128I=$((h128I+1));;
        *${myname}*NQ*128h*) h128NQ=$((h128NQ+1));;
        *${myname}*R*128i*) i128R=$((i128R+1));;
        *${myname}*I*128i*) i128I=$((i128I+1));;
        *${myname}*NQ*128i*) i128NQ=$((i128NQ+1));;
        *${myname}*R*128j*) j128R=$((j128R+1));;
        *${myname}*I*128j*) j128I=$((j128I+1));;
        *${myname}*NQ*128j*) j128NQ=$((j128NQ+1));;
        *${myname}*R*512a*) a512R=$((a512R+1));;
        *${myname}*I*512a*) a512I=$((a512I+1));;
        *${myname}*NQ*512a*) a512NQ=$((a512NQ+1));;
        *${myname}*R*512b*) b512R=$((b512R+1));;
        *${myname}*I*512b*) b512I=$((b512I+1));;
        *${myname}*NQ*512b*) b512NQ=$((b512NQ+1));;
        *${myname}*R*512c*) c512R=$((c512R+1));;
        *${myname}*I*512c*) c512I=$((c512I+1));;
        *${myname}*NQ*512c*) c512NQ=$((c512NQ+1));;
        *${myname}*R*512d*) d512R=$((d512R+1));;
        *${myname}*I*512d*) d512I=$((d512I+1));;
        *${myname}*NQ*512d*) d512NQ=$((d512NQ+1));;
        *${myname}*R*512e*) e512R=$((e512R+1));;
        *${myname}*I*512e*) e512I=$((e512I+1));;
        *${myname}*NQ*512e*) e512NQ=$((e512NQ+1));;
        *${myname}*R*2048a*) a2048R=$((a2048R+1));;
        *${myname}*I*2048a*) a2048I=$((a2048I+1));;
        *${myname}*NQ*2048a*) a2048NQ=$((a2048NQ+1));;
    esac
done < ${tmp_dir}/tmp_queue.txt

printf "\n class| R | I | NQ       class| R | I | NQ       class| R | I | NQ\n"
printf " -----+---+---+----      -----+---+---+----      -----+---+---+----\n"
printf "  32s | %d | %d | %2d       128a | %d | %d | %2d       512a | %d | %d | %2d \n" \
$s32R $s32I $s32NQ $a128R $a128I $a128NQ $a512R $a512I $a512NQ
printf "  32m | %d | %d | %2d       128b | %d | %d | %2d       512b | %d | %d | %2d \n" \
$m32R $m32I $m32NQ $b128R $b128I $b128NQ $b512R $b512I $b512NQ
printf " -----+---+---+----      128c | %d | %d | %2d       512c | %d | %d | %2d \n" \
$c128R $c128I $c128NQ $c512R $c512I $c512NQ
printf "  32a | %d | %d | %2d       128d | %d | %d | %2d       512d | %d | %d | %2d \n" \
$a32R $a32I $a32NQ $d128R $d128I $d128NQ $d512R $d512I $d512NQ
printf "  32b | %d | %d | %2d       128e | %d | %d | %2d       512e | %d | %d | %2d \n" \
$b32R $b32I $b32NQ $e128R $e128I $e128NQ $e512R $e512I $e512NQ
printf "  32c | %d | %d | %2d       128f | %d | %d | %2d       -----+---+---+----\n" \
$c32R $c32I $c32NQ $f128R $f128I $f128NQ
printf " -----+---+---+----      128g | %d | %d | %2d \n" \
$g128R $g128I $g128NQ
printf "  32d | %d | %d | %2d       128h | %d | %d | %2d \n" \
$d32R $d32I $d32NQ $h128R $h128I $h128NQ
printf "  32e | %d | %d | %2d       128i | %d | %d | %2d       -----+---+---+----\n" \
$e32R $e32I $e32NQ $i128R $i128I $i128NQ
printf "  32f | %d | %d | %2d       128j | %d | %d | %2d       2048a| %d | %d | %2d \n" \
$f32R $f32I $f32NQ $j128R $j128I $j128NQ $a2048R $a2048I $a2048NQ
printf " -----+---+---+----      -----+---+---+----      -----+---+---+----\n\n"

rm ${tmp_dir}/tmp_queue.txt