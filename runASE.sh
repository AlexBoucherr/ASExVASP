#!/bin/bash
#BATCH --account=scw1039
#SBATCH -p htc
#SBATCH -o file.ase.%J
#SBATCH -e error.%J
#SBATCH --ntasks-per-node=40      # tasks to run per node
#SBATCH -n 60                     # number of parallel processes (tasks)
#SBATCH -t 72:00:00                # time limit
#SBATCH -J TS1_vib                # Job name


# Load the environment
	module purge
	module load vasp/5.4.4-sol
	module load python/3.7.0   #python/3.5.1-chemy
	module load compiler/gnu/7/3.0
	module load ase/3.20.1

        ulimit -s unlimited
# Setting directories and defining variables
        NNODES=$SLURM_NNODES
        NCPUS=$SLURM_NTASKS
        PPN=$SLURM_NTASKS_PER_NODE
        POTPATH=$HOME/POTCAR_old/
        MYPATH=$SLURM_SUBMIT_DIR
        WDPATH=/scratch/$USER/$SLURM_JOBID
# Setting environment and creating scratch directory
        rm -rf ${WDPATH} ; mkdir -p ${WDPATH} ; cd ${WDPATH} ;
        cp ${MYPATH}/* .
        echo ${MYPATH} >> output
#**************************************************************************************
export ASE_VASP_COMMAND="mpirun vasp_std >> output"
export VASP_PP_PATH="/home/c.c1981790/VASP/PBEv54/"	# PBE from version 54
######$ export ASE_VASP_VDW=$HOME/<path-to-vdw_kernel.bindat-folder>
#**************************************************************************************
# Launch the parallel job Using ncpus processes
        env; echo VASP Start Time is `date` running NCPUs=$NCPUS PPN=$PPN
        start="$(date +%s)"
#    time mpirun vasp_std  >> output

#	which python >> output
    python ${WDPATH}/scriptVASP.py | tee vasp.out

        echo VASP Finish Time is `date` ; stop="$(date +%s)" ; finish=$(( $stop-$start ))
        echo VASP $SLURM_JOBID Job-Time  $finish seconds
# Copy output data to home
        rm CHG PCDAT XDATCAR* IBZKPT PI* REPORT
        mv ${WDPATH}/* ${MYPATH}/
        rm -rf ${WDPATH}
