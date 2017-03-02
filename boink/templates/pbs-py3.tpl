cd $PBS_O_WORKDIR

module unload python
module load GNU/4.8.3
module load anaconda

source activate py3.partitioning
