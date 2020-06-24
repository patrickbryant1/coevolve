#!/bin/bash
#SBATCH -A SNIC2019-3-319
#SBATCH -c 4
#SBATCH -t 2:00:00
#SBATCH --array=1-50
#SBATCH --error=/home/p/pbryant/pfs/results/coevolve/bench/trrosetta/yeast/pos/error/%A_%a.error
#SBATCH --output=/home/p/pbryant/pfs/results/coevolve/bench/trrosetta/yeast/pos/out/%A_%a.out


#Run trRosetta
MSADIR=/home/j/juliezhu/pfs/coevolve_yeast/data_coev/intact/pos/combined_data/msa
PROTEINPAIRS=/home/p/pbryant/pfs/runs/coevolve/trrosetta/test/yeast/yeast_pos.txt

SINGIMG=/home/j/juliezhu/pfs/trRosetta/trRosetta-ipython.simg
ROSETTA=/proj/nobackup/snic2019-35-62/TrRosetta
PAIR=$(sed -n $SLURM_ARRAY_TASK_ID'p' $PROTEINPAIRS)
echo $PAIR
BASENAME=$MSADIR/$PAIR
OUTDIR=/home/p/pbryant/pfs/results/coevolve/bench/trrosetta/yeast/pos/
#Insert flexibility in MSA
SEQLENS=/home/p/pbryant/pfs/coevolve/src/seqlen.txt
/opt/singularity/3.5.3/bin/singularity exec -B /proj:/proj,/home:/home /pfs/nobackup/home/p/pbryant/singularity/bio.sif python3 /home/p/pbryant/pfs/coevolve/src/model_rosetta/insert_flex.py --infile $BASENAME.a3m --aa_pattern 'G' --pattern_length 20 --seqlens $SEQLENS --outdir $OUTDIR

#Predict distance and torsion angles
/opt/singularity/3.5.3/bin/singularity run -B /proj:/proj,/home:/home $SINGIMG python $ROSETTA/trRosetta_1/network/predict.py -m $ROSETTA/model2019_07 $OUTDIR$PAIR'_flexed.a3m' $OUTDIR$PAIR.npz
#The npz file contains the  distances and torison angles

#Step 2 - construct a pdb model from the .npz file
#NOTE a fasta file containing the sequence is needed
#Make fasta
#head -n 2 $BASENAME.a3m > $BASENAME.fasta
#Make pdb
#/opt/singularity3/bin/singularity run --nv /home/pbryant/singularity_ims/trRosetta-gpu.simg python /home/pbryant/trRosetta/trRosetta2/trRosetta.py $BASENAME.npz  $BASENAME.fasta $BASENAME.pdb
