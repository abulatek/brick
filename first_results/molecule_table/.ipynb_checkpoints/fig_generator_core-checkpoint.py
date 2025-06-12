import os
os.chdir("/blue/adamginsburg/abulatek/brick/first_results/molecule_table/data_model_figs")
import glob
specfns_c = sorted(glob.glob("*_c_*.pdf"))
for fn in specfns_c:
    print(f"\\includegraphics[width=\\textwidth]{{figures/spectra/core/{fn}}}\\\\")
os.chdir("/blue/adamginsburg/abulatek/brick/first_results/molecule_table")