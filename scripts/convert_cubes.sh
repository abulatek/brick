#!/bin/bash
#SBATCH --job-name=convert_cubes        # Job name
#SBATCH --output convert_cubes_%j.log   # Output and error log name
#SBATCH --mail-user=abulatek@ufl.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1					    # Run on a single CPU
#SBATCH --mem=20gb					    # Job memory request
#SBATCH --time=24:00:00					# Wall time limit HH:MM:SS
#SBATCH --qos=adamginsburg-b			# Allocation

cd /blue/adamginsburg/abulatek/brick

pwd; hostname; date

echo "Running cube conversion script on 1 CPU core"

/blue/adamginsburg/abulatek/conda/envs/py312/bin/python3 /blue/adamginsburg/abulatek/brick/first_results/reduction/3_convert_cubes_w_pbcor.py

date