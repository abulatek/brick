########################################################
# Started Logging At: 2022-02-02 17:30:42
########################################################
########################################################
# # Started Logging At: 2022-02-02 17:30:43
########################################################
########################################################
# Started Logging At: 2022-02-02 17:30:50
########################################################
########################################################
# # Started Logging At: 2022-02-02 17:30:51
########################################################
########################################################
# Started Logging At: 2022-02-02 17:31:58
########################################################
########################################################
# # Started Logging At: 2022-02-02 17:31:58
########################################################
# Synthetic spectrum
get_ipython().run_line_magic('run', '/blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py')
########################################################
# Started Logging At: 2022-02-02 17:33:44
########################################################
########################################################
# # Started Logging At: 2022-02-02 17:33:45
########################################################
# Synthetic spectrum
get_ipython().run_line_magic('run', '/blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py')
########################################################
# Started Logging At: 2022-02-02 17:35:21
########################################################
########################################################
# # Started Logging At: 2022-02-02 17:35:22
########################################################
# Synthetic spectrum
get_ipython().run_line_magic('run', '/blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py')
########################################################
# Started Logging At: 2022-02-02 17:35:40
########################################################
########################################################
# # Started Logging At: 2022-02-02 17:35:41
########################################################
# Synthetic spectrum
get_ipython().run_line_magic('run', '/blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py')
# Synthetic spectrum
get_ipython().run_line_magic('run', '/blue/adamginsburg/abulatek/brick/first_results/LTE/reworked_LTE_lineforest.py')
fmin = 110330.34*u.MHz # ch3cncube.spectral_axis.min() # I should not have to hard-code these...
fmax = 110383.6*u.MHz # ch3cncube.spectral_axis.max()
ch3cn_freqs, ch3cn_A, ch3cn_g, ch3cn_E_U, ch3cn_partfunc = get_molecular_parameters('CH3CN', 
                                                                                     fmin=fmin, 
                                                                                     fmax=fmax, 
                                                                                     catalog='JPL')
# used to be 'CH3CN' as first argument
print(ch3cn_E_U)
