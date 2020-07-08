#!/usr/bin/env bash
IP=./Y2H_union.txt
NF5=./nf5.txt
OUTDIR=./
./analyze_bench.py --interacting_pairs $IP --nf5_pairs $NF5 --outdir $OUTDIR
