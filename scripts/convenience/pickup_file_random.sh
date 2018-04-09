#!/bin/bash

### N.file or N.dir in each directory should be less than 32767.

function rand_pickup() {
    ls_dir=$1
    dir_lst=(`ls $ls_dir`)
    rand_dir=${dir_lst[$(($RANDOM % ${#dir_lst[@]}))]}
    if [ -f $ls_dir/$rand_dir ]; then
	echo $ls_dir/$rand_dir
    else
	rand_pickup $ls_dir/$rand_dir
    fi
}

a=`rand_pickup .`; echo $a
