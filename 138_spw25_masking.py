## Pylab requires this to have a visual component, so if using pylab, this should be ran from an xpra GUI session
from spectral_cube import SpectralCube
from astropy import units as u

print("Now masking the cube at sigma threshold")
# Use the code below to mask the cube at a certain sigma threshold
cube = SpectralCube.read('source_ab_138_spw25_dirty_512.image')
mask = cube > 7.0*u.mJy/u.beam

from casatools import image
ia = image()
boolmask = mask.include().compute()

print("Now outputting the mask file")
# Use the code below to output a mask file
ia.open('source_ab_138_spw25_dirty_512.image')
cs = ia.coordsys()
ia.close()

ia.fromarray(outfile='source_ab_138_spw25_dirty_512.mask', pixels=boolmask.astype('float')[:,None,:,:].T, csys=cs.torecord(), overwrite=True)
ia.close()

print("Now visualizing masked cube without the circular vignette")
# Use the code below to visualize what the cube looks like without its circular vignette
import pylab as pl; pl.ion(); pl.draw(); pl.show()
# unmasked = cube.with_mask(np.ones(cube.shape, dtype='bool'), inherit_mask=False)
# Above sometimes crashes Python. Let's try: 
unmasked = cube.with_mask(cube == cube, inherit_mask=False)

# print("Now visualizing the erosion and dilation")
# Use the code below to visualize the erosion and dilation
import scipy.ndimage
# x = boolmask[2228] 
# pl.figure()
# pl.imshow(x)
# y = scipy.ndimage.binary_erosion(x, iterations=2)
# pl.imshow(y)
# z1 = scipy.ndimage.binary_dilation(scipy.ndimage.binary_erosion(x, iterations=2), iterations=2)
# pl.imshow(z1)
# z2 = scipy.ndimage.binary_dilation(scipy.ndimage.binary_erosion(x, iterations=3), iterations=4)
# pl.imshow(z2)

print("Now calculating the eroded/dilated mask")
# Use the code below to output an eroded/dilated mask
boolmask_e3_d4 = scipy.ndimage.binary_dilation(scipy.ndimage.binary_erosion(boolmask, iterations=3), iterations=4)
print("Now outputting the eroded/dilated mask")
ia.fromarray(outfile='source_ab_138_spw25_dirty_512_e3_d4.mask', pixels=boolmask_e3_d4.astype('float')[:,None,:,:].T, csys=cs.torecord(), overwrite=True)
ia.close()
