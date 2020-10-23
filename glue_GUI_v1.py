from spectral_cube import SpectralCube, StokesSpectralCube # maybe don't need StokesSpectralCube
from glue.core import Data, DataCollection
from glue.app.qt.application import GlueApplication
from glue.viewers.image.qt import ImageViewer
from glue.core.coordinates import coordinates_from_wcs

image_dir = '/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/'
image_file = image_dir + 'BrickMaser_87_spw25.residual'
cube = SpectralCube.read(image_file, format = 'casa_image')

# Try to replicate what spectral_cube_to_data(cube) does
label = None
cube = StokesSpectralCube({'I': cube})
result = Data(label=label)
result.coords = coordinates_from_wcs(cube.wcs)

for component in cube.components:
    data = getattr(cube, component)._data
    result.add_component(data, label='STOKES {0}'.format(component))

result._preferred_translation = SpectralCube

# create a GUI session
dc = DataCollection(result)
ga = GlueApplication(dc)

image = ga.new_data_viewer(ImageViewer)
image.add_data(result) # did not take long to run

# show the GUI
ga.start()

