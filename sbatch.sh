#!/bin/bash

#SBATCH -J test3           # Job name
#SBATCH -o test3.o%j       # Name of stdout output file
#SBATCH -e test3.e%j       # Name of stderr error file
#SBATCH -p spr             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for OpenMP)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for OpenMP)
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=yjung3@tacc.utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

# Other commands must follow all #SBATCH directives...

module load python

# module list
pwd
date

# Set thread count (default value is 1)...

export OMP_NUM_THREADS=48   # this is 1 thread/core; may want to start lower

# Always run your jobs out of $SCRATCH.  Your input files, output files, 
# and exectuable should be in the $SCRATCH directory hierarchy.  
# Change directories to your $SCRATCH directory where your executable is

cp test3.py $SCRATCH
cp /home1/99999/yjung3/filamentFields/filamentFields.cpython-39-x86_64-linux-gnu.so $SCRATCH
cd $SCRATCH

# Launch OpenMP code...
python test3.py