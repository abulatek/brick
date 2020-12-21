from spectral_cube import SpectralCube
import pyspeckit as sp

results_path = '/blue/adamginsburg/abulatek/brick/symlinks/imaging_results/'

cube = SpectralCube.read(results_path+'source_ab_138_spw31_clean_2sigma_n100000_masked_3sigma.image', format='casa_image')
print(cube.shape)
subcube = cube[2927:3036]

spec = sp.Spectrum(subcube)

# Note to self: this probably won't work
