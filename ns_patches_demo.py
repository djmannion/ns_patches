
import os.path
import collections
import random

import numpy as np

import psychopy.visual, psychopy.event, psychopy.filters

import stimuli.mcgill_db


def main():

	img_list = stimuli.mcgill_db.local_db_info()

	random.shuffle( img_list )

	ap = collections.namedtuple( "Aperture", [ "ecc",
	                                           "theta",
	                                           "diam",
	                                           "xy",
	                                           "type",
	                                           "hemi",
	                                           "id",
	                                           "ring"
	                                         ]
	                           )

	ring_thetas = [ [ 0, 60, 120, 180, 240, 300 ],
	                [ 0, 35, 70, 110, 145, 180, 215, 250, 290, 325 ],
	                [ 0, 30, 60, 120, 150, 180, 210, 240, 300, 330 ],
	              ]

	ring_ecc = [ 62, 123, 212 ]

	ring_diam = [ 89, 117, 158 ]

	mask_info = []

	mask_count = 0

	for ( i_ring, ( thetas, ecc, diam ) ) in enumerate( zip( ring_thetas,
	                                                         ring_ecc,
	                                                         ring_diam
	                                                       )
	                                                  ):

		for theta in thetas:

			xy = psychopy.misc.pol2cart( theta, ecc )

			if theta < 90 or theta > 270:
				hemi = "R"
			else:
				hemi = "L"

			mask_info.append( ap( ecc = ecc,
			                      theta = theta,
			                      diam = diam,
			                      xy = xy,
			                      type = "",
			                      hemi = hemi,
			                      id = mask_count,
			                      ring = i_ring
			                    )
			                )

			mask_count += 1

	masks = make_masks( mask_info )

	n_masks = len( masks )

	n_coh = 13

	win = psychopy.visual.Window( ( 512, 512 ),
	                              fullscr = False,
	                              allowGUI = True,
	                              units = "pix"
	                            )

	try:
		stims = init_stims( win, masks )

		keep_going = True

		while keep_going:

			set_stims( win, n_coh, img_list, stims )

			_ = [ stim.draw() for stim in stims ]

			win.flip()

			keys = psychopy.event.waitKeys()

			if "q" in keys:
				keep_going = False
			elif "s" in keys:
				win.getMovieFrame()
				win.saveMovieFrames( "cap.png" )


	except:
		win.close()
		raise

	win.close()


def set_stims( win, n_coh, img_list, stims ):

	coh_img = random.choice( img_list )

	other_img = [ random.choice( img_list ) for _ in stims ]

	img = [ coh_img ] * n_coh + other_img[ :( len( stims ) - n_coh ) ]
	random.shuffle( img )

	for ( i_stim, stim ) in enumerate( stims ):

		print img[ i_stim ]

		img_dir, img_name = os.path.split( img[ i_stim ] )

		img_path = os.path.join( "/home/damien/mcgill_db",
		                         img_dir, img_name
		                       )

		img_rast = stimuli.mcgill_db.read_image( img_path, True, True, True )

		stim.setImage( np.flipud( img_rast[ :512, :512 ] ) )

	return stims

def init_stims( win, masks ):

	stim = collections.namedtuple( "stim", [ "tex", "img", "i_mask", "type" ] )

	stims = []

	n_stim = len( masks )

	for mask in masks:

		stims.append( psychopy.visual.ImageStim( win = win,
		                                         image = np.ones( ( 512, 512 ) ),
		                                         mask = mask,
		                                         units = "pix",
		                                         size = ( 512, 512 ),
		                                       )
		            )

	return stims


def make_masks( mask_info ):

	masks = []

	for mask in mask_info:

		radius = mask.diam / 2.0 / 512
		xy = np.array( mask.xy ) / ( 512.0 / 2.0 )

		mask_rast = psychopy.filters.makeMask( matrixSize = 512,
		                                       shape = "raisedCosine",
		                                       radius = radius,
		                                       center = xy
		                                     )

		masks.append( mask_rast )

	return masks


if __name__ == "__main__":
	main()
