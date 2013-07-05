
import os.path
import csv

import numpy as np
import scipy.stats

import psychopy.misc

import ns_patches.paths


class ConfigContainer( object ):
	pass


def get_conf( subj_id = None ):

	conf = ConfigContainer()

	conf.stim = _get_stim_conf()
	conf.acq = _get_acq_conf()
	conf.exp = _get_exp_conf( conf )
	conf.loc = _get_loc_conf( conf )
	conf.ana = _get_ana_conf()
	conf.all_subj = _get_subj_conf()

	if subj_id is not None:
		conf.subj = _get_subj_conf( subj_id )

	return conf


def _get_exp_conf( conf ):

	exp_conf = ConfigContainer()

	exp_conf.id = "ns_patches"

	exp_conf.mod_patches = range( 28 )
	_ = [ exp_conf.mod_patches.remove( p )
	      for p in [ 20 - 1, 26 - 1 ]
	    ]

	exp_conf.n_mod_patches = len( exp_conf.mod_patches )

	exp_conf.n_incoh_patches = exp_conf.n_mod_patches / 2

	exp_conf.n_img = 20

	exp_conf.n_img_rep_per_run = 4

	exp_conf.n_img_incoh_per_run = 2

	exp_conf.n_trials = exp_conf.n_img * exp_conf.n_img_rep_per_run

	exp_conf.n_null = exp_conf.n_trials / 4

	exp_conf.n_seq = exp_conf.n_trials + exp_conf.n_null

	exp_conf.n_pre = 8

	exp_conf.bin_len_s = 4

	exp_conf.n_run_seq = exp_conf.n_seq + exp_conf.n_pre

	exp_conf.n_censor_vols = ( exp_conf.n_pre *
	                           exp_conf.bin_len_s /
	                           conf.acq.tr_s
	                         )

	exp_conf.n_vol = ( ( exp_conf.n_seq + exp_conf.n_pre ) *
	                   exp_conf.bin_len_s /
	                   conf.acq.tr_s
	                 )

	exp_conf.run_len_s = exp_conf.n_vol * conf.acq.tr_s

	exp_conf.img_on_s = 2.0

	return exp_conf


def _get_ana_conf():

	ana_conf = ConfigContainer()

	return ana_conf

def _get_acq_conf():
	"""Get the acquisition configuration"""

	acq_conf = ConfigContainer()

	acq_conf.monitor_name = "UMN_7T"
	acq_conf.test_monitor_name = "N13_CRT"

	acq_conf.tr_s = 2.0

	# how to reshape the data to be in +RAS convention
	# subsequent commands are relative to the data AFTER this operation
	# see docs for how to determine these
	acq_conf.ras = ( "-x", "-z", "-y" )

	# phase encode direction, according to the data's internal axes
	acq_conf.ph_enc_dir = "x"

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
	acq_conf.dwell_ms = 0.72 / 2.0

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
	                [ 30, 150, 210, 330 ]
	              ]

	# ~ 1.8, 3.5, 6.1, 8.5 deg
	ring_ecc = [ 62, 123, 212, 340 ]

	# M0 = 29.2/3.67; M = 29.2/(e+3.67); D = M0/M*15*2
	ring_diam = [ 45, 59, 80, 99 * 340./297 ]

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


	stim_conf.mask_edge = "raisedCos"
	stim_conf.mask_edge_frac = 0.4

	stim_conf.img_diam_pix = 1024

	return stim_conf

def _get_subj_conf( subj_id = None ):

	s1000 = ConfigContainer()

	s1000.subj_id = "s1000"
	s1000.subj_id_loc = "s1000_loc"
	s1000.acq_date = ""
	s1000.acq_date_loc = "20130620"
	s1000.comments = ""
	s1000.n_runs = 12
	s1000.n_loc_runs = 6
	s1000.mot_base = 7
	s1000.mot_base_loc = 4

	s1000.extra_al_params = [ "-parang", "1", "-13", "-3",
	                          "-parang", "2", "16", "26",
	                          "-parang", "3", "6", "16",
	                          "-maxrot", "10",
	                          "-source_automask+2",
	                          "-nocmass"
	                        ]

	s1000.extra_al_params_loc = [ "-parang", "1", "-13", "-3",
	                              "-parang", "2", "16", "26",
	                              "-parang", "3", "6", "16",
	                              "-maxrot", "10",
	                              "-source_automask+2",
	                              "-nocmass"
	                            ]

	s1000.node_k = { "lh" : 130318,
	                 "rh" : 131151
	               }


	subj = ConfigContainer()

	subj.subj = { "s1000" : s1000,
	            }

	if subj_id is None:
		return subj
	else:
		return subj.subj[ subj_id ]


def make_loc_timing( conf ):

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

		# loop over runs
		for _ in xrange( 1, conf.loc.n_max_runs + 1 ):

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


def check_loc_timing( conf, run_timing ):

	cmd = [ "3dDeconvolve",
	        "-nodata",
	        "{n:d}".format( n = conf.loc.n_vol ),
	        "{n:0f}".format( n = conf.acq.tr_s ),
	        "-polort", "A",
	        "-CENSORTR", "0-{n:d}".format( n = conf.loc.n_cull_vol - 1 ),
	        "-num_stimts", "{n:d}".format( n = conf.stim.n_patches ),
	        "-local_times",
	      ]

	for ( i_patch, patch ) in enumerate( run_timing ):

		cmd.extend( [ "-stim_times",
		              "{n:d}".format( n = i_patch + 1 ),
		              "'1D: " + " ".join( map( str, patch ) ) + "'",
		              "'SPMG1(1)'"
		            ]
		          )

	print " ".join( cmd )


def gen_exp_patch_timing( conf ):
	"""Generates the image timing patterns for each patch.

	Paramters
	---------
	conf : ns_patches.config.get_conf() object
		Configuration info

	Returns
	-------
	design : numpy array of ints, ( n_trials, n_patches )
		Each entry contains the image index for each (modulated) patch and trial

	"""

	rows = np.array( [ [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
	                   [ 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1 ],
	                   [ 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0 ],
	                   [ 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1 ]
	                 ]
	               )

	base_design = np.zeros( ( conf.exp.n_trials, conf.exp.n_mod_patches ) )

	for i_img in xrange( conf.exp.n_img ):

		patch_order = np.random.permutation( conf.exp.n_mod_patches )

		i_design = ( np.arange( conf.exp.n_img_rep_per_run  ) +
		             i_img * conf.exp.n_img_rep_per_run
		           )

		# temporarily have images be one-based for multiplication purposes
		base_design[ i_design, : ] = ( i_img + 1 ) * rows[ :, patch_order ]

	# return to zero-based image indexing
	# 'null' events are now -1
	base_design -= 1

	success = False

	na = np.repeat( np.arange( conf.exp.n_img ), conf.exp.n_img_rep_per_run / 2 )

	while not success:

		design = base_design.copy()

		extras = np.tile( np.arange( conf.exp.n_img ),
		                  conf.exp.n_img_rep_per_run / 2
		                )

		for i_patch in xrange( conf.exp.n_mod_patches ):

			i_blanks = np.where( design[ :, i_patch ] < 0 )[ 0 ]

			good = False

			while not good:

				i_extras = np.random.permutation( len( extras ) )

				fillers = extras[ i_extras ]

				good = np.all( fillers != na )

			design[ i_blanks, i_patch ] = fillers

		success = np.all( scipy.stats.mode( design, axis = 1 )[ 1 ] ==
		                  conf.exp.n_incoh_patches
		                )

	design = design[ np.random.permutation( conf.exp.n_trials ), : ]

	# all images equally represented
	assert np.all( [ np.sum( design == x ) ==
	                 conf.exp.n_mod_patches * conf.exp.n_img_rep_per_run
	                 for x in xrange( conf.exp.n_img )
	               ]
	             )

	( mode_img, mode_k ) = scipy.stats.mode( design, axis = 1 )

	# half incoherent patches per trial
	assert np.all( mode_k == conf.exp.n_incoh_patches )

	# each image is the mode the correct number of times
	assert np.all( [ np.sum( mode_img == x ) == conf.exp.n_img_rep_per_run
	                 for x in xrange( conf.exp.n_img )
	               ]
	             )

	return design
