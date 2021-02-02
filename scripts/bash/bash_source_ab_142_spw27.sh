#!/bin/bash
#SBATCH --job-name=bash_ab_142_spw27_3sigma_pbmask0p18	# Job name
#SBATCH --mail-user=abulatek@ufl.edu			# Where to send mail
#SBATCH --mail-type=END,FAIL				# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output ab_142_spw27_3sigma_pbmask0p18_%j.log	# Output and error log name
#SBATCH --nodes=1					# Maximum number of nodes to be allocated
#SBATCH --ntasks=1					# Run on a single CPU
#SBATCH --cpus-per-task=4				# Number of CPU cores per task
#SBATCH --mem=20gb					# Job memory request
#SBATCH --time=24:00:00					# Wall time limit HH:MM:SS
#SBATCH --qos=adamginsburg-b				# Allocation

pwd; hostname; date

export LOGFILENAME='casa_clean_ab_142_spw27_3sigma_pbmask0p18.log'

echo $LOGFILENAME

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa
export CASA_newpath=/orange/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa

xvfb-run -d ${CASA_newpath} --logfile=${LOGFILENAME}  --nogui --nologger -c "execfile('/blue/adamginsburg/abulatek/brick/scripts/cleans/clean_source_ab_142_spw27.py')"
# ${CASA_newpath} -c --nogui --nologger /blue/adamginsburg/abulatek/brick/clean_script.py	# This is likely not the right command

