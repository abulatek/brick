#!/bin/bash
#SBATCH --job-name=imaging_152_spw31	  # Job name
#SBATCH --mail-user=abulatek@ufl.edu	  # Where to send mail
#SBATCH --mail-type=END,FAIL			  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output imaging_152_spw31_%j.log # Output and error log name
#SBATCH --nodes=1					      # Maximum number of nodes to be allocated
#SBATCH --ntasks=1					      # Run on a single CPU
#SBATCH --cpus-per-task=4				  # Number of CPU cores per task
#SBATCH --mem=64gb					      # Job memory request
#SBATCH --time=24:00:00					  # Wall time limit HH:MM:SS
#SBATCH --qos=adamginsburg-b			  # Allocation

pwd; hostname; date

export LOGFILENAME='casa_clean_ab_152_spw31_2sigma.log'

echo $LOGFILENAME

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa
export CASA6=/blue/adamginsburg/adamginsburg/casa/casa-6.1.0-118/bin/casa
export CASA66=/blue/adamginsburg/adamginsburg/casa/casa-6.6.0-2-py3.8.el8/bin/casa

xvfb-run -d ${CASA66} --logfile=${LOGFILENAME}  --nogui --nologger -c "execfile('/blue/adamginsburg/abulatek/brick/scripts/imaging/imaging_source_ab_152_spw31.py')"

