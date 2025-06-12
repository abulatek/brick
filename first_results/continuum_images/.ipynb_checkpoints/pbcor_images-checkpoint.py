# sbatch --job-name=pbcor_59_5d_61 --output=pbcor_59_5d_61-%j.log --account=astronomy-dept --qos=astronomy-dept-b --ntasks=8 --mail-type=ALL --mail-user=abulatek@ufl.edu --nodes=1 --mem=8gb --time=4:00:00 --constraint=el8 --wrap "/blue/adamginsburg/adamginsburg/casa/casa-6.6.0-2-py3.8.el8/bin/casa -c 'execfile(\"/blue/adamginsburg/abulatek/brick/first_results/continuum_images/pbcor_images.py\")'"

impbcor(imagename="brickmaser_cont_data_59_mtmfs_incl_adjparam.image.tt0", 
        pbimage="brickmaser_cont_data_59_mtmfs_incl_adjparam.pb.tt0", 
        outfile="brickmaser_cont_data_59_mtmfs_incl_adjparam.image.tt0.pbcor")

impbcor(imagename="brickmaser_cont_data_5d_mtmfs_incl_adjparam.image.tt0", 
        pbimage="brickmaser_cont_data_5d_mtmfs_incl_adjparam.pb.tt0", 
        outfile="brickmaser_cont_data_5d_mtmfs_incl_adjparam.image.tt0.pbcor")

impbcor(imagename="brickmaser_cont_data_61_mtmfs_incl_adjparam.image.tt0", 
        pbimage="brickmaser_cont_data_61_mtmfs_incl_adjparam.pb.tt0", 
        outfile="brickmaser_cont_data_61_mtmfs_incl_adjparam.image.tt0.pbcor")