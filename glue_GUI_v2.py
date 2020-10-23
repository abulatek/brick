from spectral_cube import SpectralCube, StokesSpectralCube
from glue.core import Data, DataCollection
from glue.app.qt.application import GlueApplication
from glue.viewers.image.qt import ImageViewer
from glue.core.coordinates import coordinates_from_wcs

### Import the .image, .model, and .residual files
image_file = 'BrickMaser_87_spw25.image'
model_file = 'BrickMaser_87_spw25.model'
resid_file = 'BrickMaser_87_spw25.residual'

## Get the full file paths
file_dir = '/orange/adamginsburg/brick_alma_linesurvey/2019.1.00092.S/imaging_results/'
image_path = file_dir + image_file
model_path = file_dir + model_file
resid_path = file_dir + resid_file

## Import the files as cubes using SpectralCube
print("Import files with SpectralCube")
image_cube = SpectralCube.read(image_path, format = 'casa_image')
model_cube = SpectralCube.read(model_path, format = 'casa_image')
resid_cube = SpectralCube.read(resid_path, format = 'casa_image')

## Replicate what spectral_cube_to_data(cube) does
for cube in [image_cube, model_cube, resid_cube]: # I hope this works
    label = None
    cube = StokesSpectralCube({'I': cube})

## Initialize results
print("Initialize results")
image_result = Data(label='image')
model_result = Data(label='model')
resid_result = Data(label='resid')

## Get WCS coordinates for files
print("Get WCS coordinates")
image_result.coords = coordinates_from_wcs(image_cube.wcs)
model_result.coords = coordinates_from_wcs(model_cube.wcs)
resid_result.coords = coordinates_from_wcs(resid_cube.wcs)

## Add components, from Adam's code
print("Add cube components to results")
image_result.add_component(image_cube, 'image_cube')
model_result.add_component(model_cube, 'model_cube')
resid_result.add_component(resid_cube, 'resid_cube')

## Ensure SpectralCube is used
for result in [image_result, model_result, resid_result]:
    result._preferred_translation = SpectralCube

## Create a GUI session
print("Create a GUI session")
dc = DataCollection([image_result, model_result, resid_result])
ga = GlueApplication(dc)

print("Add data to viewer")
view = ga.new_data_viewr(ImageViewer)
view.add_data(image_result)
view.add_data(model_result)
view.add_data(resid_result)

