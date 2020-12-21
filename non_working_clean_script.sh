#!/bin/bash
#SBATCH --job-name=ab_138_spw25_2sigma_clean
#SBATCH --mail-user=abulatek@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output ab_138_spw25_2sigma_clean_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --qos=adamginsburg-b
pwd; hostname; date

echo $LOGFILENAME

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa

xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger -c "execfile('/blue/adamginsburg/abulatek/brick/clean_script.py')"

