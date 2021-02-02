#!/bin/bash
#SBATCH --job-name=testing_bash
#SBATCH --mail-user=abulatek@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output testing_bash_script_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=00:01:00

echo "This is Alyssa, testing if she can write a bash script"
