
import numpy as np

import psychopy.misc


class ConfigContainer( object ):
	pass


def get_conf( subj_id = None ):

	conf = ConfigContainer()

	conf.exp = _get_exp_conf()
	conf.loc = _get_loc_conf()
	conf.stim = _get_stim_conf()

	return conf


def _get_exp_conf():

	exp_conf = ConfigContainer()

	exp_conf.id = "ns_patches"

	return exp_conf


def _get_loc_conf():

	loc_conf = ConfigContainer()

	loc_conf.id = "ns_patches_loc"

	return loc_conf


def _get_stim_conf():

	stim_conf = ConfigContainer()

	stim_conf.n_patches = 32

	patch_dt = np.dtype( [ ( "id", "int" ),
	                       ( "ecc", "float" ),
	                       ( "theta", "float" ),
	                       ( "cx", "float" ),
	                       ( "cy", "float" ),
	                       ( "vf", "|S1" ),
	                       ( "ring", "int" )
	                     ]
	                   )

	stim_conf.patches = np.empty( ( stim_conf.n_patches ), dtype = patch_dt )

	ring_thetas = [ [ 0, 60, 120, 180, 240, 300 ],
	                [ 0, 35, 70, 110, 145, 180, 215, 250, 290, 325 ],
	                [ 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 ],
	                [ 45, 135, 225, 315 ]
	              ]

	# ~ 1.8, 3.5, 6.1 deg
	ring_ecc = [ 62, 123, 212, 262 ]

	# M0 = 29.2/3.67; M = 29.2/(e+3.67); D = M0/M*30*2
	ring_diam = [ 89, 117, 158, 183 ]

	i_patch = 0

	for ( i_ring, ( thetas, ecc, diam ) ) in enumerate( zip( ring_thetas,
	                                                         ring_ecc,
	                                                         ring_diam
	                                                       )
	                                                  ):

		for theta in thetas:

			patch = stim_conf.patches[ i_patch ]

			( cx, cy ) = psychopy.misc.pol2cart( theta, ecc )

			patch[ "id" ] = i_patch
			patch[ "ecc" ] = ecc
			patch[ "theta" ] = theta
			patch[ "cx" ] = cx
			patch[ "cy" ] = cy
			patch[ "ring" ] = i_ring

			if cx < 0:
				patch[ "vf" ] = "L"
			else:
				patch[ "vf" ] = "R"

			i_patch += 1


	return stim_conf




