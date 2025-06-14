#! /bin/bash

if [ $# -ne 1 ]; then
    echo "Error: No input file is specified."
    exit 1
fi

#python pocket_trk.py -f 12 -sig AFSD -prn 2-8 -sig AFSP -prn 2-8 -IQ $1 
python pocket_trk.py -f 12 -sig AFSD -prn 8 -IQ -p $1 

