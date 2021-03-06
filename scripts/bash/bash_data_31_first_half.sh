#!/bin/bash
#SBATCH --job-name=imaging				# Job name
#SBATCH --output imaging_%j.log                         # Output and error log name
#SBATCH --mail-user=abulatek@ufl.edu			# Where to send mail
#SBATCH --mail-type=END,FAIL				# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1					# Maximum number of nodes to be allocated
#SBATCH --ntasks=1					# Run on a single CPU
#SBATCH --cpus-per-task=4				# Number of CPU cores per task
#SBATCH --mem=40gb					# Job memory request
#SBATCH --time=24:00:00					# Wall time limit HH:MM:SS
#SBATCH --qos=adamginsburg-b				# Allocation

pwd; hostname; date

export CASA=/blue/adamginsburg/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa
export CASA6=/blue/adamginsburg/adamginsburg/casa/casa-6.1.0-118/bin/casa

export data_dir=data_31
export field='2'

freq_list=('144' '146' '127' '129' '137' '139' '137')
spw_list=('49' '51' '65' '67' '69' '71' '85')

for i in "${!freq_list[@]}"
do
    export freq_num=${freq_list[i]}
    export spw=${spw_list[i]}
    export LOGFILENAME='casa_clean_ab_'$freq_num'_spw'$spw'_2sigma.log'
    echo $data_dir 
    echo $freq_num 
    echo $spw 
    echo $field
    echo $LOGFILENAME
    ${CASA6} --logfile=${LOGFILENAME} --nogui --nologger -c '/blue/adamginsburg/abulatek/brick/scripts/imaging/imaging.py' $data_dir $freq_num $spw $field
#    xvfb-run -d ${CASA6} --logfile=${LOGFILENAME}  --nogui --nologger -c "execfile('/blue/adamginsburg/abulatek/brick/scripts/imaging/imaging.py')"
done

