#!/bin/bash

delete_d=("__pycache__")
delete_f=(".DS_Store" "*~" ".*~" "._*" "*.pyc")

for d in ${delete_d[@]}; do
    find . -type d -name $d -print -exec rm -fr {} \;
done

for f in ${delete_f[@]}; do
    find . -type f -name $f -print -exec rm -f  {} \;
done
