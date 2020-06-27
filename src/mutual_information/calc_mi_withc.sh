INFILE=O13297_Q01159_100.a3m
OUTDIR=/home/patrick/results/coevolve/test/C++/
#complile
g++ -Wall -std=c++11 -O cmi.cpp -o fastMI
wait
#run
./fastMI $INFILE $OUTDIR
