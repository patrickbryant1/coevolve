#!/bin/bash
#SBATCH -A SNIC2019-3-319
#SBATCH -c 2
#SBATCH -t 10:00:00
#SBATCH --array=1-800
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

## load modules for singularity
ml GCC/8.2.0-2.31.1  OpenMPI/3.1.3 ADIOS2/2.5.0-Python-3.7.2
module load singularity/2.6.1

abspath=/home/j/juliezhu/pfs/coevolve_yeast
alipath=/proj/nobackup/snic2019-35-62/juliezhu
nflst=$abspath/nf5.txt


##SLURM_ARRAY_TASK_ID=2

## find the zst file name: 
beginID=$(python3 -c "print((${offset}+(${SLURM_ARRAY_TASK_ID}-1))*671+1)")
endID=$(python3 -c "print((${offset}+${SLURM_ARRAY_TASK_ID})*671)")
echo $beginID,$endID

##make tmpdir to unzip the zst file
tmpdir=$(mktemp -d)

cd $tmpdir
cp $alipath/nf5/${beginID}_${endID}.tar.zst ./

/home/p/pbryant/pfs/zstd/programs/zstd -d ${beginID}_${endID}.tar.zst
tar -xvf ${beginID}_${endID}.tar

#mktmp 
tmpdir=$(mktemp -d)


for INFILE in `ls *.a3m`;do 
##run cal_mi_withc.sh 
echo infile is $INFILE
#wc -l $INFILE 
pref=$(echo $INFILE | cut -d. -f1)
#info=$(/home/g/gabriele/pfs/coevolve/src/mutual_information/fastMI $INFILE $alipath/mi_nf5/ | grep $pref )
singularity exec --nv $abspath/gdca_tf.sif python3 $abspath/src/compute_score.py $INFILE $tmpdir/$pref >> $alipath/miscore_nf5.txt
#/opt/singularity/2.6.1/bin/singularity exec --nv $abspath/gdca_tf.sif $abspath/src/compute_score.py $INFILE $alipath/mi_nf5/$pref >> $abspath/miscore_nf5.txt
done

cd $tmpdir
tar --use-compress-program /home/p/pbryant/pfs/zstd/programs/zstd -cf ${SLURM_ARRAY_TASK_ID}.tar.zst ./*
mv ${SLURM_ARRAY_TASK_ID}.tar.zst $alipath/mi_nf5/
##clean tmpdir
rm -rf $tmpdir



