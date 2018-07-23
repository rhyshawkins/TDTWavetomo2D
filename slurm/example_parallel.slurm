#!/bin/bash
#SBATCH --job-name=example_parallel
#SBATCH --partition=transcale
#SBATCH --time=1:00:00
###########################################################
# USER PARAMETERS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-socket=14
###########################################################
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000mb
#SBATCH -o output_%j.txt
#SBATCH -e stderr_%j.txt
#
# Intel 17.0.4
source /softs/intel/l_ics/2017_update4/compilers_and_libraries_2017.4.196/linux/bin/compilervars.sh intel64
source /softs/intel/l_ics/2017_update4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpivars.sh intel64
export INTEL_LICENSE_FILE=/softs/intel/l_ics/license:$INTEL_LICENSE_FILE

#
# Binary Location 
#

BASE=/home_nfs/lgltpe/rhys.hawkins/software/TDTWavetomo2D/
BIN=$BASE/wavetomo2dfrequencysliceinvert_pt

#
# Input observations
#
DATABASE=/erc_transcale/rhys.hawkins/data/
INPUT=$DATABASE/data20.txt
PRIOR=$DATABASE/prior.txt

#
# Output directory
#
OUTPUTBASE=/erc_transcale/rhys.hawkins/wavetomo2d
mkdir -p $OUTPUTBASE

OUTPUT=$OUTPUTBASE/example_parallel/
mkdir -p $OUTPUT

#
# Command line parameters
#
OPTIONS="-i $INPUT \
	-o $OUTPUT \
        -M $PRIOR \
	-x 6 -y 6 -w 4 \
	-s 0 \
        -n 64 -N 95.5 \
        -a 19 -A 50.5 \
        -k 500 \
        -t 10000 \
        -v 100 \
	-E \
        -c 7 \
	" 


# To run an mpi code with Omnipath network
mpirun -PSM2 -n $SLURM_NTASKS -ppn $SLURM_CPUS_ON_NODE $BIN $OPTIONS

exit
