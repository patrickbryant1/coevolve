#!/bin/bash
#SBATCH -A SNIC2019-3-319
#SBATCH -c 4
#SBATCH -t 96:00:00
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

#load modules
ml GCC/8.2.0-2.31.1  OpenMPI/3.1.3 ADIOS2/2.5.0-Python-3.7.2
module load singularity/2.6.1

##get the argument from the input
if [ -z $1 ]
then
        offset=0
else
        offset=$1
fi

##abspath
abspath=/home/j/juliezhu/pfs/coevolve_yeast

#interacting pairs list
#posfile=/home/j/juliezhu/pfs/coevolve_yeast/ppi_yeast/uniprot_ppi.txt
posfile=/home/j/juliezhu/pfs/coevolve_yeast/ppi_yeast/ppi_neg.txt

idx=$(expr $offset + $SLURM_ARRAY_TASK_ID)
#idx=17
#get the two proteins
pr1=$(sed -n ${idx}p $posfile | cut -f1)
pr2=$(sed -n ${idx}p $posfile | cut -f2)
echo $pr1 $pr2

##get the homolog sequences 
## blast orthologs
file1=$abspath/ortholog/$pr1
file2=$abspath/ortholog/$pr2

## make sure that two homolog groups exist
if [[ ! -s $file1 || ! -s $file2 ]]; then
echo "one or more is missing!"
exit 1
fi

#make tmpdir
tmpdir=$(mktemp -d)
#tmpdir=/scratch/jtest
echo $tmpdir

mkdir $tmpdir/seq_uniprot
mkdir $tmpdir/tmp
mkdir $tmpdir/homolog_seq
mkdir $tmpdir/jackhmmer
mkdir $tmpdir/combined_data
mkdir $tmpdir/combined_data/seq
mkdir $tmpdir/combined_data/msa

## make sure the pair do not interact
#if grep -q "${pr1}.*${pr2}\|${pr2}.*${pr1}" $posfile;then
#echo 'they are interacting pair!'
#exit 1
#fi

##extract fasta file 
grep -A 1 $pr1 $abspath/data_raw/yeast/UP000002311_559292.fasta > $tmpdir/seq_uniprot/$pr1.fa
grep -A 1 $pr2 $abspath/data_raw/yeast/UP000002311_559292.fasta > $tmpdir/seq_uniprot/$pr2.fa

##run python code to keep species intersection
##generate two-column species info at the same time
python3 $abspath/src/intersection.py $file1 $file2 $idx $tmpdir/tmp
ls $tmpdir/tmp
#run blastdbcmd to retrieve sequences for homologs
cd $tmpdir/homolog_seq/
blastdbcmd -db /home/j/juliezhu/pfs/coevolve/data_raw/refome/blast_db/allproteome -entry_batch $tmpdir/tmp/${idx}_${pr1} -out ${idx}_${pr1}.fa
blastdbcmd -db /home/j/juliezhu/pfs/coevolve/data_raw/refome/blast_db/allproteome -entry_batch $tmpdir/tmp/${idx}_${pr2} -out ${idx}_${pr2}.fa
echo "blastdbcmd is done!"

##one liner for homolog file.
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr1}.fa )" > ${idx}_${pr1}.fa
sed -i '1d' ${idx}_${pr1}.fa
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr2}.fa )" > ${idx}_${pr2}.fa
sed -i '1d' ${idx}_${pr2}.fa

## run jackhmmer
cd $tmpdir/jackhmmer/
jackhmmer -N 1 -A ${idx}_${pr1}.sto $tmpdir/seq_uniprot/$pr1.fa $tmpdir/homolog_seq/${idx}_${pr1}.fa 
jackhmmer -N 1 -A ${idx}_${pr2}.sto $tmpdir/seq_uniprot/$pr2.fa $tmpdir/homolog_seq/${idx}_${pr2}.fa 
##reformat from sto to a3m
/home/j/juliezhu/pfs/programs/hhsuite2/scripts/reformat.pl ${idx}_${pr1}.sto ${idx}_${pr1}.a3m
/home/j/juliezhu/pfs/programs/hhsuite2/scripts/reformat.pl ${idx}_${pr2}.sto ${idx}_${pr2}.a3m
#ls $tmpdir/jackhmmer
rm ${idx}_${pr1}.sto
rm ${idx}_${pr2}.sto

##one liner for msa in a3m format
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr1}.a3m )" > ${idx}_${pr1}.a3m
sed -i '1d' ${idx}_${pr1}.a3m
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr2}.a3m )" > ${idx}_${pr2}.a3m
sed -i '1d' ${idx}_${pr2}.a3m

##trim the a3m to have the same length
/home/j/jlamb/pfs/call_scripts/a3mToTrimmed.py ${idx}_${pr1}.a3m > ${idx}_${pr1}.a3m_trim
/home/j/jlamb/pfs/call_scripts/a3mToTrimmed.py ${idx}_${pr2}.a3m > ${idx}_${pr2}.a3m_trim
#ls $tmpdir/jackhmmer
mv ${idx}_${pr1}.a3m_trim ${idx}_${pr1}.a3m
mv ${idx}_${pr2}.a3m_trim ${idx}_${pr2}.a3m
echo "jackhmmer MSA(whole process) is done!"

## trim the output of jackhmmer. find the missing sequences in the jackhmmer process and trim them both in their original homolog fasta files and output msa files. so that the fasta files and msa files could match to each other. 
## seq_path is the path to homolog sequences.
## msa_path is the path to output of jackhmmer.
## eco_path is the path to ecoli sequences.

## reorder a3m
python $abspath/src/reorder_a3m.py $tmpdir/homolog_seq/${idx}_${pr1}.fa $tmpdir/homolog_seq/${idx}_${pr2}.fa $tmpdir/jackhmmer/${idx}_${pr1}.a3m $tmpdir/jackhmmer/${idx}_${pr2}.a3m

paste -d '' $tmpdir/jackhmmer/${idx}_${pr1}.a3m $tmpdir/jackhmmer/${idx}_${pr2}.a3m > $tmpdir/combined_data/msa/${pr1}_${pr2}.a3m
paste -d '' $tmpdir/seq_uniprot/$pr1.fa $tmpdir/seq_uniprot/$pr2.fa > $tmpdir/combined_data/seq/${pr1}_${pr2}.fa


cd $tmpdir/jackhmmer/
tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf ${idx}_jackhmmer.tar.zst ./*
mv ${idx}_jackhmmer.tar.zst $abspath/data_coev/intact/neg/jackhmmer/
cd $tmpdir/combined_data/
tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf ${idx}_final.tar.zst ./*
mv ${idx}_final.tar.zst $abspath/data_coev/intact/neg/combined_data/
