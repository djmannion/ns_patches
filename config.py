
import os.path
import csv

import numpy as np

import psychopy.misc

import ns_patches.paths


class ConfigContainer( object ):
	pass


def get_conf( subj_id = None ):

	conf = ConfigContainer()

	conf.exp = _get_exp_conf()
	conf.stim = _get_stim_conf()
	conf.acq = _get_acq_conf()
	conf.loc = _get_loc_conf( conf )

	return conf


def _get_exp_conf():

	exp_conf = ConfigContainer()

	exp_conf.id = "ns_patches"

	return exp_conf


def _get_acq_conf():
	"""Get the acquisition configuration"""

	acq_conf = ConfigContainer()

	acq_conf.monitor_name = "UMN_7T"

	acq_conf.tr_s = 2.0

	# how to reshape the data to be in +RAS convention
	# subsequent commands are relative to the data AFTER this operation
	# see docs for how to determine these
	acq_conf.ras = ( "-x", "-z", "-y" )

	# phase encode direction, according to the data's internal axes
	acq_conf.ph_encode_dir = "z"

	# axis index that corresponds to the inplanes
	# this is zero-based; 0 = LR, 1 = PA, 2 = IS (assuming reshape_to_RAS has been
	# set correctly)
	acq_conf.slice_axis = 1

	# direction in which slices were acquired along the axis
	acq_conf.slice_acq_dir = "+1"

	# number of slices acquired
	acq_conf.n_slices = 36

	# TE difference in the fieldmaps
	acq_conf.delta_te_ms = 1.02

	# corresponds to the echo spacing
	acq_conf.dwell_ms = 0.65 / 2.0

	return acq_conf


def _get_loc_conf( conf ):

	loc_conf = ConfigContainer()

	loc_conf.id = "ns_patches_loc"

	loc_conf.dur_s = 1.0
	loc_conf.reversal_interval_s = 0.15

	loc_conf.n_vol = 166
	loc_conf.n_cull_vol = 16
	loc_conf.n_valid_vol = loc_conf.n_vol - loc_conf.n_cull_vol

	loc_conf.run_len_s = loc_conf.n_vol * conf.acq.tr_s

	loc_conf.bin_dur_s = 4.0

	loc_conf.n_max_runs = 6

	return loc_conf


def _get_stim_conf():

	stim_conf = ConfigContainer()

	stim_conf.n_patches = 32

	patch_dt = np.dtype( [ ( "id", "int" ),
	                       ( "ecc", "float" ),
	                       ( "theta", "float" ),
	                       ( "diam", "float" ),
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

	# ~ 1.8, 3.5, 6.1, 8.5 deg
	ring_ecc = [ 62, 123, 212, 297 ]

	# M0 = 29.2/3.67; M = 29.2/(e+3.67); D = M0/M*15*2
	ring_diam = [ 45, 59, 80, 99 ]

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
			patch[ "diam" ] = diam
			patch[ "cx" ] = cx
			patch[ "cy" ] = cy
			patch[ "ring" ] = i_ring

			if cx < 0:
				patch[ "vf" ] = "L"
			else:
				patch[ "vf" ] = "R"

			i_patch += 1


	return stim_conf


def make_timing( conf ):

	exp_paths = ns_patches.paths.get_exp_paths( conf )

	n_evt = ( conf.loc.n_valid_vol *
	          conf.acq.tr_s /
	          conf.loc.bin_dur_s
	        )

	assert n_evt.is_integer()

	n_evt = int( n_evt )

	n_null_evt = n_evt / 4

	n_pre_evt = ( conf.loc.n_cull_vol *
	              conf.acq.tr_s /
	              conf.loc.bin_dur_s
	            )

	assert n_pre_evt.is_integer()

	n_pre_evt = int( n_pre_evt )

	# each patch gets its own timing
	for i_patch in xrange( conf.stim.n_patches ):

		timing_path = os.path.join( exp_paths.timing_dir,
		                            ( "ns_patches-loc_timing_patch_" +
		                              "{n:02d}.txt".format( n = i_patch )
		                            )
		                          )

		timing_file = open( timing_path, "w" )

		timing_csv = csv.writer( timing_file, delimiter = " " )

		for run_num in xrange( 1, conf.loc.n_max_runs + 1 ):

			t = np.ones( ( n_evt ) )

			t[ :n_null_evt ] = 0

			np.random.shuffle( t )

			assert t.sum() == ( n_evt - n_null_evt )

			t = np.concatenate( ( t[ -n_pre_evt: ], t ) )

			assert len( t ) == ( n_pre_evt + n_evt )
			assert ( len( t ) * conf.loc.bin_dur_s == conf.loc.run_len_s )

			t_sec = np.where( t > 0 )[ 0 ] * conf.loc.bin_dur_s

			# randomly add a TR offset. this assumes that the bin duration is double
			# the TR
			if np.random.rand() > 0.5:
				t_sec += conf.acq.tr_s

			timing_csv.writerow( [ "{n:.0f}".format( n = n ) for n in t_sec ] )

		timing_file.close()


def check_timing( conf, run_timing ):

	cmd = [ "3dDeconvolve",
	        "-nodata",
	        "{n:d}".format( conf.loc.n_vol ),
	        "{n:0f}".format( conf.acq.tr_s ),
	        "-polort", "A",
	        "-CENSORTR", "0-{n:d}".format( n = conf.loc.n_cull_vol - 1 ),
	        "-num_stimts", "{n:d}".format( n = conf.stim.n_patches ),
	        "-local_times",
	      ]

	for ( i_patch, patch ) in run_timing:

		cmd.extend( [ "-stim_times",
		              "{n:d}".format( n = i_patch + 1 ),
		              "'1D: " + " ".join( map( str, patch ) ) + "'",
		              "'SPMG1(1)'"
		            ]
		          )

	print " ".join( cmd )




