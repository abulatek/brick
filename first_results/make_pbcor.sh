#!/bin/bash
#SBATCH --job-name=make_pbcor            # Job name
#SBATCH --output make_pbcor_%j.log       # Output and error log name
#SBATCH --mail-user=abulatek@ufl.edu       # Where to send mail
#SBATCH --mail-type=END,FAIL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                          # Maximum number of nodes to be allocated
#SBATCH --ntasks=8                         # Run on a single CPU
#SBATCH --cpus-per-task=1                  # Number of CPU cores per task
#SBATCH --mem=32gb                         # Job memory request
#SBATCH --time=96:00:00                    # Wall time limit HH:MM:SS
#SBATCH --qos=adamginsburg-b               # Allocation

pwd; hostname; date

freq_spw_list=('87_spw25' '89_spw27' '91_spw25' '93_spw27' '95_spw25' '97_spw27' '98_spw29' '99_spw31' '102_spw23' '102_spw29' '104_spw25' '103_spw31' '106_spw29' '107_spw31' '110_spw29' '111_spw31' '112_spw27' '114_spw29' '127_spw65' '129_spw67' '130_spw105' '132_spw107' '134_spw45' '135_spw47' '137_spw85' '137_spw69' '139_spw71' '140_spw109' '142_spw111' '144_spw49' '146_spw51' '147_spw89' '149_spw91' '142_spw27' '152_spw31' '244_spw65' '245_spw67' '247_spw105' '249_spw107' '250_spw25' '252_spw27' '254_spw85' '255_spw87' '257_spw45' '259_spw47' '259_spw71' '261_spw109' '263_spw111' '264_spw29' '266_spw31' '268_spw89' '270_spw91' '271_spw49' '273_spw51') # This is not a totally complete list, some potential duplicates have been removed

for i in "${!freq_spw_list[@]}"
do
    export freq_spw=${freq_spw_list[i]}
    echo $freq_spw
    /blue/adamginsburg/abulatek/anaconda/bin/python3 make_pbcor.py $freq_spw False
done