#!/usr/bin/env bash
IP=./Y2H_union.txt
NF5=./nf5.txt
ID_CONV=./gene_id_meta_uniprot.tab
OUTDIR=./
./analyze_bench.py --interacting_pairs $IP --nf5_pairs $NF5 --id_conversion $ID_CONV --outdir $OUTDIR
