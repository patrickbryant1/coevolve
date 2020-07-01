#!/bin/bash
#SBATCH -A SNIC2019-3-319
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH --array=1-663
##SBATCH -p largemem
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

### collect paired alignments that have nfscore greater than 5.

abspath=/home/j/juliezhu/pfs/coevolve_yeast
alipath=/proj/nobackup/snic2019-35-62/juliezhu
comlst=$abspath/data_raw/yeast/yeast_combination.txt 
nflst=$abspath/nf5.txt 

##each parallel job computes 671 times, and there are 2663 parallel jobs which I send to array.
##similar to paired_alignment.sh
##SLURM_ARRAY_TASK_ID=3

#make tmpdir: save paired alignments.
savedir=$(mktemp -d)

####test part
##beginID=$(python3 -c "print((${offset}+(${SLURM_ARRAY_TASK_ID}-1))*10+1)")
##endID=$(python3 -c "print((${offset}+${SLURM_ARRAY_TASK_ID})*10)")

beginID=$(python3 -c "print((${offset}+(${SLURM_ARRAY_TASK_ID}-1))*671+1)")
endID=$(python3 -c "print((${offset}+${SLURM_ARRAY_TASK_ID})*671)")

echo "beginID is ${beginID},endID is ${endID},start loop!"
for i in `seq ${beginID} ${endID}`;do
echo "row in nf5.txt is ${i}"
pr1=$(sed -n ${i}p $nflst | cut -d' ' -f1)
pr2=$(sed -n ${i}p $nflst | cut -d' ' -f2)
echo "pr1 and pr2 is ${pr1},${pr2}"

#index is the row number of the pair in yeast_combination.txt, also the name when I saved the paired alignment.
index=$(grep -n "$pr1.*$pr2" $comlst | cut -d':' -f1 )

###found which folder (pw_filter${endnum}) this alignment saved in.
endnum=$(python3 -c "import math;print(math.ceil(${index}/5495000))")

fname1=$(python3 -c "print(${index}-${index}%5495+1)")
fname2=$(python3 -c "print(${index}-${index}%5495+5495)")

echo "index is ${index}, pw_filter folder is ${endnum},filename is ${fname1}_${fname2}"
##mktmp dir: unzip 5495 files and get the alignment we want. then save to $savedir.
tmpdir=$(mktemp -d)
cp $alipath/pw_filter${endnum}/${fname1}_${fname2}.tar.zst $tmpdir/

##unzip
cd $tmpdir/
/home/p/pbryant/pfs/zstd/programs/zstd -d ${fname1}_${fname2}.tar.zst
tar -xvf ${fname1}_${fname2}.tar

##save the right alignment to $savedir
cp $tmpdir/combined_msa/$index $savedir/${pr1}_${pr2}.a3m
rm -rf $tmpdir
done

cd $savedir/

tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf ${beginID}_${endID}.tar.zst ./*
mv ${beginID}_${endID}.tar.zst $alipath/nf5/
echo "job done!"
