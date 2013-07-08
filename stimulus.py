import numpy as np

import psychopy.visual, psychopy.filters

import stimuli.utils, stimuli.mcgill_db


def load_img( img_path ):

	image = stimuli.mcgill_db.read_image( img_path = img_path,
	                                      linearise = False,
	                                      convert_to_grey = False,
	                                      normalise = True
	                                    )

	# pad to a power of two, as required by psychopy (well, the mask anyway)
	image = stimuli.utils.pad_image( image )

	image = np.flipud( image )

	return image


def make_mask( patch,
               img_size,
               mask_shape = "raisedCosine",
               mask_fringe = 0.4,
               mask_range = ( -1, +1 )
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
	                  fringeWidth = mask_fringe,
	                  range = mask_range
	                )

	return mask


def get_fixation( win, conf ):

	fix_stim = []

	h_frac = 0.625

	l_s = [ ( +1, +0 ), ( +h_frac, +1 ), ( +0, +1 ), ( -h_frac, +1 ) ]
	l_e = [ ( -1, +0 ), ( -h_frac, -1 ), ( +0, -1 ), ( +h_frac, -1 ) ]

	grid = [ psychopy.visual.Line( win,
	                               start = ls,
	                               end = le,
	                               units = "norm",
	                               lineColor = [ -0.25 ] * 3,
	                               lineWidth = 1.5
	                             )
	         for ( ls, le ) in zip( l_s, l_e )
	       ]

	fix_stim.extend( grid )

	grids_r_va = [ 0.5, 1.8, 3.5, 6.1, 8.5 ]

	grid_lum = -0.25

	# list of circle stimuli
	rings = [ psychopy.visual.Circle( win,
	                                  radius = grid_r_va,
	                                  units = "deg",
	                                  edges = 96,
	                                  lineColor = grid_lum,
	                                  lineWidth = 1.5
	                                )
	          for grid_r_va in grids_r_va
	        ]

	fix_stim.extend( rings )

	outlines = [ psychopy.visual.Circle( win = win,
	                                     pos = ( patch[ "cx" ], patch[ "cy" ] ),
	                                     units = "pix",
	                                     radius = patch[ "diam" ] / 2.0,
	                                     lineColor = [ -0.25 ] * 3,
	                                     fillColor = [ 0 ] * 3,
	                                     lineWidth = 1.5,
	                                   )
	             for patch in conf.stim.patches
	           ]

	fix_stim.extend( outlines )

	fixation = psychopy.visual.Circle( win = win,
	                                   radius = 4,
	                                   units = "pix",
	                                   fillColor = [ -1, 1, -1 ],
	                                   lineWidth = 1
	                                 )

	fix_stim.append( fixation )

	return fix_stim
