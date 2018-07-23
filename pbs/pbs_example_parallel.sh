#!/bin/bash
#PBS -P xe8
#PBS -l wd
#PBS -q normal
#PBS -l walltime=7:00:00,mem=64GB,ncpus=256

HOMEDIR=/home/158/rph158
OUTPUTBASE=/short/xe8/rph158/wavetomo2d

export LONMIN=-24.50
export LONMAX=-13.50
export LATMIN=63.50
export LATMAX=66.50

export PHDBASE=$HOMEDIR/PhD/PhDS
export WAVETOMO2DBASE=$PHDBASE/tools/Wavetomo2dFrequencyInvert
export WAVETOMO2DINVERT_PT=$WAVETOMO2DBASE/wavetomo2dfrequencysliceinvert_pt
export WAVETOMO2D_OBS=$WAVETOMO2DBASE/Observations/all_joint_zz.obs
export PRIORFILE=$WAVETOMO2DBASE/raijin/laplace2dprior.txt
export LAMBDASTD=0.05
export BETASTD=0.05

export RESULTS_DIR=$OUTPUTBASE/iceland655_zz_${SLICE}_0/

export WAVETOMO2DARGS="-i $WAVETOMO2D_OBS -o $RESULTS_DIR -M $PRIORFILE -t 1000000 -x 6 -y 5 -z 5 -v 10000 -T 1 -c 16 -w 4 -k 2000 -n $LONMIN -N $LONMAX -a $LATMIN -A $LATMAX -s $SLICE -l 1.0 -H $LAMBDASTD -L $BETASTD"

mkdir -p $OUTPUTBASE
mkdir -p $RESULTS_DIR

echo $WAVETOMO2DARGS > $RESULTS_DIR/run0.args

mpirun -np $PBS_NCPUS $WAVETOMO2DINVERT_PT $WAVETOMO2DARGS 

exit
