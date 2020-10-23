import numpy as np

from glue.core import Data, DataCollection
from glue.app.qt.application import GlueApplication
try:
    from glue.viewers.image.qt.data_viewer import ImageViewer
except ImportError:
    from glue.viewers.image.qt.viewer_widget import ImageWidget as ImageViewer

from glue.core.coordinates import coordinates_from_wcs

from spectral_cube import SpectralCube


cube1 = SpectralCube.read('/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/BrickMaser_102_spw23.image', format='casa_image')
cube2 = SpectralCube.read('/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/BrickMaser_102_spw23.model', format='casa_image')
cube3 = SpectralCube.read('/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/BrickMaser_102_spw23.residual', format='casa_image')

# create some data
d1 = Data(label='image')
d2 = Data(label='model')
d3 = Data(label='residual')
d1.coords = coordinates_from_wcs(cube1.wcs)
d2.coords = coordinates_from_wcs(cube2.wcs)
d3.coords = coordinates_from_wcs(cube3.wcs)

print(f"Adding data: {cube1}")
d1.add_component(cube1, 'imcube')
print(f"Adding data: {cube2}")
d2.add_component(cube2, 'modcube')
print(f"Adding data: {cube3}")
d3.add_component(cube3, 'rescube')

print("assembling datacollection")
dc = DataCollection([d1,d2,d3])

print("create a GUI session")
ga = GlueApplication(dc)

print("imshow")
imshow = ga.new_data_viewer(ImageViewer, data=d1)
print("add d1")
imshow.add_data(d1)
print("add d2")
imshow.add_data(d2)
print("add d3")
imshow.add_data(d3)


print("show the GUI")
ga.start()
