#!/bin/sh

function calc() {
    awk "BEGIN {print $1}"
}

function make_mean_err() {
    data=($@)
    Mean=0.0
    Mean2=0.0
    for i in `seq 0 1 $((${#data[@]}-1))`; do
        Mean=`calc "$Mean  + ${data[$i]}"`
        Mean2=`calc "$Mean2 + ${data[$i]}*${data[$i]}"`
    done
    Mean=`calc "$Mean  / ${#data[@]}"`
    Mean2=`calc "$Mean2 / ${#data[@]}"`
    Error=`calc "sqrt(($Mean2 - $Mean*$Mean) / ($((${#data[@]}-1))))"`

    echo "$Mean +/- $Error"
}

