
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
	conf.ana = _get_ana_conf()
	conf.all_subj = _get_subj_conf()

	if subj_id is not None:
		conf.subj = _get_subj_conf( subj_id )

	return conf


def _get_exp_conf():

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

	return exp_conf


def _get_ana_conf():

	ana_conf = ConfigContainer()

	return ana_conf

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


	return stim_conf

def _get_subj_conf( subj_id = None ):

	s1000_loc_pilot = ConfigContainer()

	s1000_loc_pilot.subj_id = "s1000_loc_pilot"
	s1000_loc_pilot.acq_date = "20130620"
	s1000_loc_pilot.comments = ""
	s1000_loc_pilot.n_runs = 6
	s1000_loc_pilot.n_exp_runs = 0
	s1000_loc_pilot.n_loc_runs = 6
	s1000_loc_pilot.exp_runs = [ ]
	s1000_loc_pilot.loc_runs = range( 1, 7 )
	s1000_loc_pilot.n_fmaps = 1
	s1000_loc_pilot.mot_base = 4

	s1000_loc_pilot.extra_al_params = [ "-parang", "1", "-13", "-3",
	                                    "-parang", "2", "16", "26",
	                                    "-parang", "3", "6", "16",
	                                    "-maxrot", "10",
	                                    "-source_automask+2",
	                                    "-nocmass"
	                                  ]


	s1000_loc_pilot.node_k = { "lh" : 130318,
	                           "rh" : 131151
	                         }

	s1000_loc_pilot_2runs = ConfigContainer()

	s1000_loc_pilot_2runs.subj_id = "s1000_loc_pilot_2runs"
	s1000_loc_pilot_2runs.acq_date = "20130620"
	s1000_loc_pilot_2runs.comments = ""
	s1000_loc_pilot_2runs.n_runs = 2
	s1000_loc_pilot_2runs.n_exp_runs = 0
	s1000_loc_pilot_2runs.n_loc_runs = 2
	s1000_loc_pilot_2runs.exp_runs = [ ]
	s1000_loc_pilot_2runs.loc_runs = range( 1, 3 )
	s1000_loc_pilot_2runs.n_fmaps = 1
	s1000_loc_pilot_2runs.mot_base = 2

	s1000_loc_pilot_2runs.extra_al_params = [ "-parang", "1", "-13", "-3",
	                                          "-parang", "2", "16", "26",
	                                          "-parang", "3", "6", "16",
	                                          "-maxrot", "10",
	                                          "-source_automask+2",
	                                          "-nocmass"
	                                        ]


	s1000_loc_pilot_2runs.node_k = { "lh" : 130318,
	                                 "rh" : 131151
	                               }

	s1000_loc_pilot_3runs = ConfigContainer()

	s1000_loc_pilot_3runs.subj_id = "s1000_loc_pilot_3runs"
	s1000_loc_pilot_3runs.acq_date = "20130620"
	s1000_loc_pilot_3runs.comments = ""
	s1000_loc_pilot_3runs.n_runs = 3
	s1000_loc_pilot_3runs.n_exp_runs = 0
	s1000_loc_pilot_3runs.n_loc_runs = 3
	s1000_loc_pilot_3runs.exp_runs = [ ]
	s1000_loc_pilot_3runs.loc_runs = range( 1, 4 )
	s1000_loc_pilot_3runs.n_fmaps = 1
	s1000_loc_pilot_3runs.mot_base = 3

	s1000_loc_pilot_3runs.extra_al_params = [ "-parang", "1", "-13", "-3",
	                                          "-parang", "2", "16", "26",
	                                          "-parang", "3", "6", "16",
	                                          "-maxrot", "10",
	                                          "-source_automask+2",
	                                          "-nocmass"
	                                        ]


	s1000_loc_pilot_3runs.node_k = { "lh" : 130318,
	                                 "rh" : 131151
	                               }

	s1000_loc_pilot_4runs = ConfigContainer()

	s1000_loc_pilot_4runs.subj_id = "s1000_loc_pilot_4runs"
	s1000_loc_pilot_4runs.acq_date = "20130620"
	s1000_loc_pilot_4runs.comments = ""
	s1000_loc_pilot_4runs.n_runs = 4
	s1000_loc_pilot_4runs.n_exp_runs = 0
	s1000_loc_pilot_4runs.n_loc_runs = 4
	s1000_loc_pilot_4runs.exp_runs = [ ]
	s1000_loc_pilot_4runs.loc_runs = range( 1, 5 )
	s1000_loc_pilot_4runs.n_fmaps = 1
	s1000_loc_pilot_4runs.mot_base = 4

	s1000_loc_pilot_4runs.extra_al_params = [ "-parang", "1", "-13", "-3",
	                                          "-parang", "2", "16", "26",
	                                          "-parang", "3", "6", "16",
	                                          "-maxrot", "10",
	                                          "-source_automask+2",
	                                          "-nocmass"
	                                        ]


	s1000_loc_pilot_4runs.node_k = { "lh" : 130318,
	                                 "rh" : 131151
	                               }

	s1000_loc_pilot_5runs = ConfigContainer()

	s1000_loc_pilot_5runs.subj_id = "s1000_loc_pilot_5runs"
	s1000_loc_pilot_5runs.acq_date = "20130620"
	s1000_loc_pilot_5runs.comments = ""
	s1000_loc_pilot_5runs.n_runs = 5
	s1000_loc_pilot_5runs.n_exp_runs = 0
	s1000_loc_pilot_5runs.n_loc_runs = 5
	s1000_loc_pilot_5runs.exp_runs = [ ]
	s1000_loc_pilot_5runs.loc_runs = range( 1, 6 )
	s1000_loc_pilot_5runs.n_fmaps = 1
	s1000_loc_pilot_5runs.mot_base = 4

	s1000_loc_pilot_5runs.extra_al_params = [ "-parang", "1", "-13", "-3",
	                                          "-parang", "2", "16", "26",
	                                          "-parang", "3", "6", "16",
	                                          "-maxrot", "10",
	                                          "-source_automask+2",
	                                          "-nocmass"
	                                        ]


	s1000_loc_pilot_5runs.node_k = { "lh" : 130318,
	                                 "rh" : 131151
	                               }


	subj = ConfigContainer()

	subj.subj = { "s1000_loc_pilot" : s1000_loc_pilot,
	              "s1000_loc_pilot_2runs" : s1000_loc_pilot_2runs,
	              "s1000_loc_pilot_3runs" : s1000_loc_pilot_3runs,
	              "s1000_loc_pilot_4runs" : s1000_loc_pilot_4runs,
	              "s1000_loc_pilot_5runs" : s1000_loc_pilot_5runs,
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


def make_exp_timing( conf ):

	img_seq = np.repeat( np.arange( conf.exp.n_img ),
	                     conf.exp.n_img_rep_per_run
	                   )

	np.random.shuffle( img_seq )

	trials = np.tile( img_seq, ( conf.stim.n_patches, 1 ) )

	assert trials.shape == ( conf.stim.n_patches,
	                         conf.exp.n_img * conf.exp.n_img_rep_per_run
	                       )


	trial_coh = sel_coh( conf.exp.n_mod_patches,
	                     conf.exp.n_incoh_patches,
	                     conf.exp.n_trials
	                   )

	for i_incoh_patch in xrange( conf.exp.n_mod_patches ):

		i_patch = conf.exp.mod_patches[ i_incoh_patch ]

		# find the trials where this patch is designated as incoherent
		( _, incoh_trials ) = np.where( trial_coh == i_incoh_patch )

		success = False

		while not success:

			# get the image sequence for this patch and its incoherent trials
			base_seq = trials[ i_patch, incoh_trials ]

			# generate a shuffled image sequence for the incoherent trials
			perm_seq = base_seq[ np.random.permutation( len( base_seq ) ) ]

			# test whether the shuffled sequence has replaced an image with the same
			# image - no good
			diff_img = ( base_seq - perm_seq ) != 0

			if np.all( diff_img ):

				success = True

				trials[ i_patch, incoh_trials ] = perm_seq



	return ( trials, trial_coh )


def sel_coh( n_patch, n_coh, n_trials ):

	max_coh = n_trials / 2

	success = False

	while not success:

		restart = False

		coh_k = np.zeros( n_patch )

		coh = np.empty( ( n_coh, n_trials ) )

		for i_trial in xrange( n_trials ):

			p = ( 1 - ( coh_k / float( max_coh ) ) )
			p /= p.sum()

			try:
				sel = np.random.choice( n_patch,
				                        n_coh,
				                        False,
				                        p
				                      )
			except ValueError:
				restart = True
				break

			coh[ :, i_trial ] = sel

			coh_k[ sel ] += 1

		if not restart:
			success = True

	assert np.all( np.array( [ ( coh == p ).sum() for p in xrange( n_patch ) ] )
	               == max_coh
	             )

	assert np.all( np.histogram( coh.flatten(), np.arange( n_patch + 1 ) )[ 0 ]
	               == max_coh
	             )

	assert np.all( np.array( [ len( np.unique( c ) ) for c in coh.T ] )
	               == n_coh
	             )

	return coh
