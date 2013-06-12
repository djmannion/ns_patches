
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
	                [ 0, 30, 60, 120, 150, 180, 210, 240, 300, 330 ],
	                [ 0, 30, 60, 120, 150, 180, 210, 240, 300, 330 ],
	              ]

	ring_ecc = [ 53, 106, 177 ]

	ring_diam = [ 56, 72, 94 ]

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

	return masks



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
