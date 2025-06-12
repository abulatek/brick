#!/bin/bash
#SBATCH --job-name=contsub_cubes        # Job name
#SBATCH --output contsub_cubes_%j.log   # Output and error log name
#SBATCH --mail-user=abulatek@ufl.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1					    # Run on a single CPU
#SBATCH --mem=32gb					    # Job memory request
#SBATCH --time=24:00:00					# Wall time limit HH:MM:SS
#SBATCH --qos=adamginsburg-b			# Allocation

cd /blue/adamginsburg/abulatek/brick

pwd; hostname; date

echo "Running cube contsub script on 1 CPU core"

/blue/adamginsburg/abulatek/conda/envs/py312/bin/python3 /blue/adamginsburg/abulatek/brick/first_results/reduction/4_contsub_cubes.py

date