#!/bin/bash
#SBATCH -A SNIC2019-3-319
#SBATCH -c 4
#SBATCH -t 40:00:00
#SBATCH --array=1-39
####SBATCH -p largemem
#SBATCH --error=/home/j/juliezhu/pfs/coevolve_yeast/error/%A_%a.error
#SBATCH --output=/home/j/juliezhu/pfs/coevolve_yeast/out/%A_%a.out

#load modules
ml GCC/8.2.0-2.31.1  OpenMPI/3.1.3 ADIOS2/2.5.0-Python-3.7.2

##get the argument from the input
if [ -z $1 ]
then
        offset=0
else
        offset=$1
fi

abspath=/home/j/juliezhu/pfs/coevolve
dbpath=/home/j/juliezhu/pfs/coevolve/data_raw/blastdb

#idx=1
idx=$(expr $offset + $SLURM_ARRAY_TASK_ID)
proteome_list=/home/j/juliezhu/pfs/coevolve_yeast/data_raw/ID_proteomename_taxinfo.txt
pname=$(sed -n ${idx}p $proteome_list | cut -d ' ' -f1)
folder=$(sed -n ${idx}p $proteome_list | cut -d ' ' -f2)
yeast=/home/j/juliezhu/pfs/coevolve_yeast/data_raw/yeast/out_cdhit.fasta
ecolidb=/home/j/juliezhu/pfs/coevolve_yeast/data_raw/yeast/saccharomyces
echo $pname

#make tmpdir
tmpdir=$(mktemp -d)
#tmpdir=/scratch/jtest
echo $tmpdir
##make dir for reciprocal best hits: forward,forward_trimmed,backward,backward_trimmed,tmp
mkdir $tmpdir/forward
mkdir $tmpdir/forward_trimmed
mkdir $tmpdir/backward
mkdir $tmpdir/backward_trimmed
mkdir $tmpdir/tmp
mkdir $tmpdir/rbh_hit
mkdir $tmpdir/rbh_hit/manual
mkdir $tmpdir/ortholog

#-------------------------------------------------------------------------------------------------------------------
### 1. reciprocal best hit process
 
## psiblast forward

cd $tmpdir/forward
psiblast -query $yeast -db $dbpath/$pname -out $pname -num_iterations 3 -evalue 0.01 -outfmt "6 qseqid sseqid qlen slen length qcovs pident qstart qend sstart send evalue"

## remove blank lines and lines with "Search has converged!" OBS: i use '-i' means change inplace which I could only run once.
sed -i '/^Search.*/d;/^$/d' $pname

echo 'psiblast forward done!'

##run the forward filtering

python $abspath/src/blast_mainfilter.py $pname $tmpdir/forward_trimmed

echo 'psiblast forward filtering done!'

## psiblast backward

qid=$(cat $tmpdir/forward_trimmed/$pname | cut -d '	' -f2 |cut -d '|' -f2 |sort -u)

##get the fasta file of the proteome
proteome_seq=$abspath/data_raw/refome/$folder/$pname.fasta

cd $tmpdir/backward

##iteration: blast all 'target' proteins in the file

while IFS= read -r line
do

##find the sequence of the target protein $line
pat1=".*${line}.*"
pat2="^>.*"
#echo $pat1,$pat2
sequence=$(sed "0,/${pat1}/d;/${pat2}/,\$d" $proteome_seq)

##find the info line of target protein $line
id=$(grep $line $proteome_seq)
#echo $id
#echo $sequence
##save all sequences of hits in one proteome to a fa file.
echo -e "${id}\n${sequence}" >> $tmpdir/tmp/$pname.fa
done <<< "$qid"

##psiblast back
psiblast -query $tmpdir/tmp/$pname.fa -db $ecolidb  -out $pname -num_iterations 3 -evalue 0.01 -outfmt "6 qseqid sseqid qlen slen length qcovs pident qstart qend sstart send evalue"

## remove blank lines and lines with "Search has converged!" OBS: i use '-i' means change inplace which I could only run once.
sed -i '/^Search.*/d;/^$/d' $pname

#echo 'psiblast backward done!'

## run for the backward filtering
python $abspath/src/blast_mainfilter.py $pname $tmpdir/backward_trimmed

#echo 'psiblast backward filtering done!'

cd $tmp
## intersection of the forward hits and backward hits; save in $tmpdir/rbh_hit
python $abspath/src/intersec_fwbw.py $pname $tmpdir

#echo 'intersection done!'
## filtering so that one ecoli protein has one hit in each proteome
python $abspath/src/parsing_blast.py $pname $tmpdir

#echo 'one hit/ecoli done!'
#----------------------------------------------------------------------------------------------

#package rbh_hit folder back 
cd $tmpdir/rbh_hit/
tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf $pname.tar.zst ./*
mv $pname.tar.zst /home/j/juliezhu/pfs/coevolve_yeast/rbh_hit/

#----------------------------------------------------------------------------------------------
##group orthologs for each protein in the yeast proteome
ortpath=/home/j/juliezhu/pfs/coevolve_yeast/ortholog

cd $tmpdir/ortholog
for m in `seq 1 5810`;do
filename=$(sed -n ${m}p /home/j/juliezhu/pfs/coevolve_yeast/data_raw/yeast/ID_saccharomuces.txt)
touch $filename
done

#length=5
length=$(wc -l < $tmpdir/rbh_hit/$pname)

for i in `seq 1 $length`
do
yeast_row=$(sed -n ${i}p $tmpdir/rbh_hit/$pname)
yeast_name=$(sed -n ${i}p $tmpdir/rbh_hit/$pname | cut -f1 | cut -d '|' -f2)
echo $yeast_name
echo $yeast_row
#echo $eco_row,$eco_name
#grep $yeast_row $tmpdir/rbh_hit/$pname > $tmpfile2
#cat $tmpfile2
#sed -i "s/$/    $pname/" $tmpfile2
echo "${yeast_row}	${pname}" >> $tmpdir/ortholog/$yeast_name
done

for i in `find . -mindepth 1 -not -empty`;do
echo $i
yeast_name=$(echo $i | cut -d '/' -f2)
cat $i >> /home/j/juliezhu/pfs/coevolve_yeast/ortholog/$yeast_name
done
