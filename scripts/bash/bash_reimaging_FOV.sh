#!/bin/bash
#SBATCH --job-name=data_reimaging_FOV           # Job name
#SBATCH --output data_reimaging_FOV_%j.log      # Output and error log name
#SBATCH --mail-user=abulatek@ufl.edu        # Where to send mail
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                           # Maximum number of nodes to be allocated
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --cpus-per-task=4                   # Number of CPU cores per task
#SBATCH --mem=40gb                          # Job memory request
#SBATCH --time=24:00:00	                    # Wall time limit HH:MM:SS
#SBATCH --qos=adamginsburg-b				# Allocation

pwd; hostname; date

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa
export CASA6=/blue/adamginsburg/adamginsburg/casa/casa-6.1.0-118/bin/casa

data_dir_list=(data_41)
# data_dir_list=(data_55 data_5d data_31 data_41) data_31
field_list=('3')
# field_list=('2' '2' '2' '3')

# freq_spw_05_04 = '93_spw27'
# freq_spw_06_05 = '110_spw29'
# freq_spw_07_06 = '130_spw105'
# freq_spw_08_07 = '146_spw51'
# freq_spw_14_13 = '257_spw45', which is actually wrong
# freq_spw_14_13 = '259_spw47'

freq_list=('259')
# freq_list=('93' '110' '146' '257') '130'
spw_list=('47')
# spw_list=('27' '29' '51' '45') '105'

for i in "${!freq_list[@]}"
do
    export data_dir=${data_dir_list[i]}
    export field=${field_list[i]}
    export freq_num=${freq_list[i]}
    export spw=${spw_list[i]}
    export LOGFILENAME='casa_clean_ab_'$freq_num'_spw'$spw'_2sigma_reimaging_FOV.log'
    echo $data_dir 
    echo $freq_num 
    echo $spw 
    echo $field
    echo $LOGFILENAME
    ${CASA6} --logfile=${LOGFILENAME} --nogui --nologger -c '/blue/adamginsburg/abulatek/brick/scripts/imaging/reimaging_FOV.py' $data_dir $freq_num $spw $field
#    xvfb-run -d ${CASA6} --logfile=${LOGFILENAME}  --nogui --nologger -c "execfile('/blue/adamginsburg/abulatek/brick/scripts/imaging/reimaging_FOV.py')"
done
