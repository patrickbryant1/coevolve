#!/bin/bash
#SBATCH -A SNIC2019-3-319
#SBATCH -c 4
#SBATCH -t 120:00:00
#SBATCH --array=1-1000
###SBATCH -p largemem
#SBATCH --error=/home/j/juliezhu/pfs/coevolve_yeast/error/%A_%a.error
#SBATCH --output=/home/j/juliezhu/pfs/coevolve_yeast/out/%A_%a.out

##get the argument from the input
if [ -z $1 ]
then
        offset=0
else
        offset=$1
fi

#load modules for python
ml GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
ml Python/3.6.6

#load modules for hmmer
ml icc/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
ml ifort/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
ml HMMER/3.2.1

##abs path and file path
abspath=/home/j/juliezhu/pfs/coevolve_yeast
file=$abspath/data_raw/yeast/yeast_combination.txt

#writeFile=/proj/nobackup/snic2019-35-62/juliezhu/ID_nf_pairaln.txt
#writeFile=$abspath/ID_nf_pairaln.txt
#writeFile=$abspath/ID_nf_pairaln2.txt
writeFile=$abspath/ID_nf_pairaln3.txt

#SLURM_ARRAY_TASK_ID=3

#make tmpdir
tmpdir=$(mktemp -d)
#echo $tmpdir

mkdir $tmpdir/seq_uniprot
mkdir $tmpdir/tmp
mkdir $tmpdir/homolog_seq
mkdir $tmpdir/jackhmmer
mkdir $tmpdir/combined_msa

##each parallel job computes 157*35=5495 times, and there are 83*37=3071 parallel jobs which I send to array.
##I need a base number so that each base number*offset equals to Nth 5495.

beginID=$(python3 -c "print((${offset}+(${SLURM_ARRAY_TASK_ID}-1))*5495+1)")
endID=$(python3 -c "print((${offset}+${SLURM_ARRAY_TASK_ID})*5495)")

##test part
#beginID=$(python3 -c "print((${offset}+(${SLURM_ARRAY_TASK_ID}-1))*10+1)")
#endID=$(python3 -c "print((${offset}+${SLURM_ARRAY_TASK_ID})*10)")

#echo "beginID is ${beginID},endID is ${endID},start loop!"
##loop through all 5495 files 
for idx in `seq ${beginID} ${endID}`;do
echo "idx is ${idx}"
#idx=$(expr $offset + $SLURM_ARRAY_TASK_ID)
#for idx in `seq 1 61`;do
line=$(sed -n ${idx}p $file)
pr1=$(echo $line | cut -d ' ' -f1)
pr2=$(echo $line | cut -d ' ' -f2)
echo $pr1,$pr2

## blast orthologs
file1=$abspath/ortholog/$pr1
file2=$abspath/ortholog/$pr2
#echo $file1,$file2

## exit current iteration if any ortholog file is empty
if [[ ! -s $file1 || ! -s $file2 ]]; then
#echo "one or more is missing!"
echo "${pr1} ${pr2} No 0 0-0 0 0 0" >> $writeFile
continue
fi


##intersection based on species
python3 $abspath/src/intersection.py $file1 $file2 $idx $tmpdir/tmp
numintersect=$(wc -l < $tmpdir/tmp/${idx}_${pr1})
##exit if there's no intersection
if [ $(wc -l < ${tmpdir}/tmp/${idx}_${pr1}) -eq 0 ];then
#echo "no intersectant sequences from ${pr1} and ${pr2}! "
echo "$pr1 $pr2 Yes 0 0-0 0 0 0" >> $writeFile
continue
fi

#echo "${pr1} and ${pr2} continue"

##grep yeast sequence
grep -A 1 $pr1 $abspath/data_raw/yeast/UP000002311_559292.fasta > $tmpdir/seq_uniprot/$pr1.fa
grep -A 1 $pr2 $abspath/data_raw/yeast/UP000002311_559292.fasta > $tmpdir/seq_uniprot/$pr2.fa

cd $tmpdir/homolog_seq
blastdbcmd -db /home/j/juliezhu/pfs/coevolve/data_raw/refome/blast_db/allproteome -entry_batch $tmpdir/tmp/${idx}_${pr1} -out ${idx}_${pr1}.fa
blastdbcmd -db /home/j/juliezhu/pfs/coevolve/data_raw/refome/blast_db/allproteome -entry_batch $tmpdir/tmp/${idx}_${pr2} -out ${idx}_${pr2}.fa


##one liner for homolog file.
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr1}.fa )" > ${idx}_${pr1}.fa
sed -i '1d' ${idx}_${pr1}.fa
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr1}.fa )" > ${idx}_${pr1}.fa
sed -i '1d' ${idx}_${pr1}.fa

## run jackhmmer
cd $tmpdir/jackhmmer/
jackhmmer -N 1 --noali --cpu 4 -A ${idx}_${pr1}.sto $tmpdir/seq_uniprot/$pr1.fa $tmpdir/homolog_seq/${idx}_${pr1}.fa
jackhmmer -N 1 --noali --cpu 4 -A ${idx}_${pr2}.sto $tmpdir/seq_uniprot/$pr2.fa $tmpdir/homolog_seq/${idx}_${pr2}.fa

/home/j/juliezhu/pfs/programs/hhsuite2/scripts/reformat.pl ${idx}_${pr1}.sto ${idx}_${pr1}.a3m
/home/j/juliezhu/pfs/programs/hhsuite2/scripts/reformat.pl ${idx}_${pr2}.sto ${idx}_${pr2}.a3m
rm ${idx}_${pr1}.sto
rm ${idx}_${pr2}.sto

##hhfilter remove highly similar sequences: pairwised alignments (Baker)
/home/j/juliezhu/pfs/programs/hhsuite2/bin/hhfilter -v 0 -id 90 -cov 75 -i ${idx}_${pr1}.a3m -o ${idx}_${pr1}.a3m
/home/j/juliezhu/pfs/programs/hhsuite2/bin/hhfilter -v 0 -id 90 -cov 75 -i ${idx}_${pr2}.a3m -o ${idx}_${pr2}.a3m

numHHfilter1A=$(wc -l < ${idx}_${pr1}.a3m)
numHHfilter1B=$(wc -l < ${idx}_${pr2}.a3m)

##one line
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr1}.a3m )" > ${idx}_${pr1}.a3m
sed -i '1d' ${idx}_${pr1}.a3m
echo "$(awk '/^>/{print a;a="";print;next}{a=a$0}END{print a}' ${idx}_${pr2}.a3m )" > ${idx}_${pr2}.a3m
sed -i '1d' ${idx}_${pr2}.a3m

#trim the a3m to have the same length
/home/j/jlamb/pfs/call_scripts/a3mToTrimmed.py ${idx}_${pr1}.a3m > ${idx}_${pr1}.a3m_trim
mv ${idx}_${pr1}.a3m_trim ${idx}_${pr1}.a3m
/home/j/jlamb/pfs/call_scripts/a3mToTrimmed.py ${idx}_${pr2}.a3m > ${idx}_${pr2}.a3m_trim
mv ${idx}_${pr2}.a3m_trim ${idx}_${pr2}.a3m

##reorder a3m so that sequences from same speices can match
python $abspath/src/reorder_a3m.py $tmpdir/homolog_seq/${idx}_${pr1}.fa $tmpdir/homolog_seq/${idx}_${pr2}.fa $tmpdir/jackhmmer/${idx}_${pr1}.a3m $tmpdir/jackhmmer/${idx}_${pr2}.a3m

numreorder=$(wc -l < $tmpdir/jackhmmer/${idx}_${pr1}.a3m)
#echo "reorder is done!"

## check if any a3m is empty after previous filtering
if [[ ! -s ${tmpdir}/jackhmmer/${idx}_${pr1}.a3m || ! -s ${tmpdir}/jackhmmer/${idx}_${pr2}.a3m ]];then
#echo "${pr1} and ${pr2} have no common msa after the filtering"
echo"$pr1 $pr2 Yes $numintersect ${numHHfilter1A}-${numHHfilter1B} 0 0 0" >> $writeFile
continue
fi

###concatenate alignments 
paste -d '' $tmpdir/jackhmmer/${idx}_${pr1}.a3m $tmpdir/jackhmmer/${idx}_${pr2}.a3m > $tmpdir/combined_msa/$idx
#
### filter the informative sequences
/home/j/juliezhu/pfs/programs/hhsuite2/bin/hhfilter -v 0 -id 90 -i $tmpdir/combined_msa/${idx} -o $tmpdir/combined_msa/${idx}
N90=$(wc -l < $tmpdir/combined_msa/${idx})
L=$(sed -n 2p $tmpdir/combined_msa/${idx} | wc -c)
L=$(bc <<< "scale=3; sqrt(($L-1))")
nf=$(python3 -c "print(round(${N90}/${L},3))")

echo "$pr1 $pr2 Yes $numintersect ${numHHfilter1A}-${numHHfilter1B} $numreorder $N90 $nf" >> $writeFile

##clear the output file if no error
numerror=$(find $abspath/error -name "*_${SLURM_ARRAY_TASK_ID}.error" | xargs grep -iF 'error' | wc -l)
if [ $numerror -eq 0 ];then
#echo "No error!"
find $abspath/out -name "*_${SLURM_ARRAY_TASK_ID}.out" | xargs truncate -s 0
fi

done
cd $tmpdir/
tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf ${beginID}_${endID}.tar.zst combined_msa/
mv ${beginID}_${endID}.tar.zst /proj/nobackup/snic2019-35-62/juliezhu/pw_filter3/
#else
#        echo "one or more is missing"
#        echo "$pr1 $pr2 0" >> /proj/nobackup/snic2019-35-62/juliezhu/ID_nf_pairaln.txt
#fi
 
