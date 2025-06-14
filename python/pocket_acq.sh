#! /bin/bash

if [ $# -ne 1 ]; then
    echo "Error: No input file is specified."
    exit 1
fi

#python pocket_acq.py -sig AFSD -prn 1-8 -f 12 -fi 0 $1 
python pocket_acq.py -sig AFSD -prn 8 -f 12 -fi 0 -p $1 

