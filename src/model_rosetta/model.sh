#!/bin/bash -l

#Run trRosetta to model protein complexes
#Step 1 - predict distances and torison angles

MSAS=/home/pbryant/coevolve/trRosetta_modeling/ecoli/pos/msa
PROTEINPAIRS=/home/pbryant/coevolve/trRosetta_modeling/ecoli_pos_lst
for i in {1..39}
do
	PAIR=$(sed -n $i'p' $PROTEINPAIRS)	
	echo $PAIR
	BASENAME=$MSAS/$PAIR
	/opt/singularity3/bin/singularity run --nv /home/pbryant/singularity_ims/trRosetta-gpu.simg python /home/pbryant/trRosetta/network/predict.py -m /home/pbryant/trRosetta/model2019_07 $BASENAME.a3m $BASENAME.npz
	#The npz file contains the  distances and torison angles

	#Step 2 - construct a pdb model from the .npz file
	#NOTE a fasta file containing the sequence is needed
	#Make fasta
	head -n 2 $BASENAME.a3m > $BASENAME.fasta
	#Make pdb
	/opt/singularity3/bin/singularity run --nv /home/pbryant/singularity_ims/trRosetta-gpu.simg python /home/pbryant/trRosetta/trRosetta2/trRosetta.py $BASENAME.npz  $BASENAME.fasta $BASENAME.pdb
done


