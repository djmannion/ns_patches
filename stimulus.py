import numpy as np

import psychopy.visual, psychopy.filters

import stimuli.utils, stimuli.mcgill_db


class Stimulus( object ):
	"""Brightness matching experiment stimulus object.

	Parameters
	----------
	win : PsychoPy window instance
		Window to draw to.
	img_file : string
		Path to the image.
	patch_info : ns_patches.config.stim.patches structure
		Info on each patch.
	active_patches : list of integers
		ID for each active patch for this image, optional.

	"""

	def __init__( self,
	              win,
	              img_file,
	              patch_info,
	              mask_edge = "raisedCosine",
	              mask_edge_frac = 0.4,
	              active_patches = None
	            ):

		self._win = win

		self._img_file = img_file

		self._patch_info = patch_info
		self._active_patches = active_patches

		self._mask_edge = mask_edge
		self._mask_edge_frac = mask_edge_frac

		if self._active_patches is None:
			self._active_patches = []

		self._set_image()

		self._set_mask()


	def draw( self ):
		"""Draws the stimulus to the window"""

		self._img_tex.draw()


	def _set_image( self ):
		"""Sets the current texture image"""

		image = stimuli.mcgill_db.read_image( img_path = self._img_file,
		                                      linearise = True,
		                                      normalise = True
		                                    )

		# pad to a power of two, as required by psychopy (well, the mask anyway)
		image = stimuli.utils.pad_image( image )

		image = np.flipud( image )

		self._img_size = image.shape

		self._img_tex = psychopy.visual.ImageStim( win = self._win,
		                                           image = image,
		                                           size = self._img_size,
		                                           units = "pix"
		                                         )

	def reset_curr_mask( self ):

		self._img_tex.setMask( self._mask )


	def disable_mask( self ):

		self._img_tex.setMask( np.ones( self._mask.shape ) )


	def enable_mask( self ):

		self._img_tex.setMask( self._mask )


	def _set_mask( self ):
		"""Sets the visibility mask"""

		active_patches = [ patch for patch in self._patch_info
		                   if patch[ "id" ] in self._active_patches
		                 ]

		# start off with the mask in [ 0, 1 ] range, to make it easier
		mask = np.zeros( ( self._img_size ) ) + 0.0

		# stop the lines below getting a bit unwieldy
		mask_func = psychopy.filters.makeMask

		for patch in active_patches:

			norm_rad = patch[ "diam" ] / self._img_size[ 0 ]
			norm_xy = ( np.array( [ patch[ "cx" ], patch[ "cy" ] ] ) /
			            ( self._img_size[ 0 ] / 2.0 )
			          )

			mask += mask_func( matrixSize = self._img_size[ 0 ],
			                   shape = self._mask_edge,
			                   radius = [ norm_rad ] * 2,
			                   center = norm_xy,
			                   fringeWidth = self._mask_edge_frac,
			                   range = [ 0, 1 ]
			                 )

		# convert to [ -1, +1 ]
		self._mask = mask * 2.0 - 1.0

		self._mask = np.clip( self._mask, -1, +1 )

		self._img_tex.setMask( self._mask )
