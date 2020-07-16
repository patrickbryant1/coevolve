#!/bin/bash
#SBATCH -A SNIC2020-5-300
#SBATCH -c 28
#SBATCH -t 10:00:00
#SBATCH --array=1-1000
##SBATCH -p largemem
#SBATCH --error=/home/j/juliezhu/pfs/coevolve_yeast/error/%A_%a.error
#SBATCH --output=/home/j/juliezhu/pfs/coevolve_yeast/out/%A_%a.out
#load modules for python
ml GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
ml Python/3.6.6

#load modules for hmmer
ml icc/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
ml ifort/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
ml HMMER/3.2.1

##get the argument from the input
if [ -z $1 ]
then
        offset=0
else
        offset=$1
fi

idx=$(expr $offset + $SLURM_ARRAY_TASK_ID)
yeastID=/home/j/juliezhu/pfs/coevolve_yeast/data_raw/yeast/ID_saccharomuces.txt
yeastproteome=/home/j/juliezhu/pfs/coevolve_yeast/data_raw/yeast/out_cdhit.fasta 
allproteome=/home/j/juliezhu/pfs/coevolve/data_raw/allproteome/whole.fasta


##mk tmpdir: save jackhmmer MSA 
tmpdir=$(mktemp -d)

## extract corresponding yeast sequence 
yname=$(sed -n ${idx}p $yeastID)
grep -A 1 $yname $yeastproteome > $tmpdir/$yname.fa

cd $tmpdir
jackhmmer -N 3 -E 0.01 --cpu 28 -A $yname.sto $tmpdir/$yname.fa $allproteome 

##a3m format
/home/j/juliezhu/pfs/programs/hhsuite2/scripts/reformat.pl $yname.sto $yname.a3m
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${yname}.a3m )" > $yname.a3m
sed -i '1d' $yname.a3m
/home/j/jlamb/pfs/call_scripts/a3mToTrimmed.py $yname.a3m > ${yname}_trimmed
mv ${yname}_trimmed $yname.a3m

##record how many sequences in the MSA
numseq=$(grep -c '>' ${yname}.a3m )
echo $yname $numseq >> /home/j/juliezhu/pfs/coevolve_yeast/result/jackhmmer_direct/jhdirect_numhit.txt
#rm fasta sequence
rm $tmpdir/$yname.fa
tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf $yname.tar.zst ./*
mv $yname.tar.zst /home/j/juliezhu/pfs/coevolve_yeast/jackhmmer_directhit/







