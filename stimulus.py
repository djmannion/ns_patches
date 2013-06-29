import numpy as np

import psychopy.visual, psychopy.filters

import stimuli.utils, stimuli.mcgill_db


class Stimulus( object ):
	"""NS patches experiment stimulus object.

	Parameters
	----------
	win : PsychoPy window instance
		Window to draw to.
	img : string or ndarray
		Either a path to an image, or the image itself.
	mask : ndarray, optional
		Mask to apply to the image.

	"""

	def __init__( self,
	              win,
	              img,
	              mask = None
	            ):

		self._win = win
		self._img = img
		self._mask = mask

		self._gen_tex()

		self.update_image( self._img )

		self.update_mask( self._mask )


	def draw( self ):
		"""Draws the stimulus to the window"""

		self._img_tex.draw()


	def _gen_tex( self ):
		"""Generates a texture image"""

		self._img_tex = psychopy.visual.ImageStim( win = self._win,
		                                           image = np.ones( ( 256, 256 ) ),
		                                           size = ( 256, 256 ),
		                                           units = "pix"
		                                         )

	def update_image( self, img ):

		self._img = img

		if type( self._img ) == np.ndarray:
			image = self._img
		else:
			image = load_img( self._img )

		self._img_size = image.shape

		self._img_tex.setImage( image )
		self._img_tex.setSize( self._img_size )


	def update_mask( self, mask ):

		self._mask = mask

		if self._mask is None:
			self._mask = np.ones( self._img_size )

		self._img_tex.setMask( mask )


def load_img( img_path ):

	image = stimuli.mcgill_db.read_image( img_path = img_path,
	                                      linearise = True,
	                                      normalise = True
	                                    )

	# pad to a power of two, as required by psychopy (well, the mask anyway)
	image = stimuli.utils.pad_image( image )

	image = np.flipud( image )

	return image


def make_mask( patch,
               img_size,
               mask_shape = "raisedCosine",
               mask_fringe = 0.4
             ):

	mask_func = psychopy.filters.makeMask

	norm_rad = patch[ "diam" ] / img_size

	norm_xy = ( np.array( [ patch[ "cx" ], patch[ "cy" ] ] ) /
	            ( img_size / 2.0 )
	          )

	mask = mask_func( matrixSize = img_size,
	                  shape = mask_shape,
	                  radius = [ norm_rad ] * 2,
	                  center = norm_xy,
	                  fringeWidth = mask_fringe
	                )

	return mask






