#!/bin/bash
#SBATCH --job-name=first_clean
#SBATCH --mail-user=abulatek@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output first_clean_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=00:30:00
pwd; hostname; date

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa

xvfb-run -d ${CASA} --logfile=${LOGFILENAME}  --nogui --nologger -c "execfile('/blue/adamginsburg/abulatek/brick/first_clean_script.py')"
