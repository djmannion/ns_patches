
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

	exp_conf.n_run_seq = exp_conf.n_seq + exp_conf.n_pre

	exp_conf.bin_len_s = 4

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

		i_design = np.arange( 4 ) + i_img * 4

		# temporarily have images be one-based for multiplication purposes
		base_design[ i_design, : ] = ( i_img + 1 ) * rows[ :, patch_order ]

	# return to zero-based image indexing
	# 'null' events are now -1
	base_design -= 1

	success = False

	na = np.repeat( np.arange( conf.exp.n_img ), 2 )

	while not success:

		design = base_design.copy()

		extras = np.tile( np.arange( conf.exp.n_img ), 2 )

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


def gen_exp_patch_timing( conf ):

	rows = np.zeros( ( 4, conf.exp.n_mod_patches )  )
	rows[ 0, :conf.exp.n_incoh_patches ] = 1
	rows[ 1, conf.exp.n_incoh_patches: ] = 1
	rows[ 2, np.arange( 0, conf.exp.n_mod_patches, 2 ) ] = 1
	rows[ 3, np.arange( 1, conf.exp.n_mod_patches, 2 ) ] = 1

	design = np.zeros( ( conf.exp.n_trials, conf.exp.n_mod_patches ) )

	for i_img in xrange( conf.exp.n_img ):

		patch_order = np.random.permutation( conf.exp.n_mod_patches )

		i_design = np.arange( 4 ) + i_img * 4

		design_unfulfilled[ i_design, : ] = i_img *  rows[ :, patch_order ]

	oops = True

	while oops:

		design = design_unfulfilled.copy()

		extras = np.tile( np.arange( 20 )[ :, np.newaxis ], ( 1, 2 ) )

		for i_patch in xrange( conf.exp.n_mod_patches ):

			design[ blanks, i_patch ] = extras[ np.random.permutation( len( extras ) ) ]

		for i_img in xrange( 20 ):
			for i_r in xrange( 4 ):
				row = i_img * 4 + i_r

	return rows


def gen_exp_patch_timing_old( conf ):
	"""Generates the image timing patterns for each patch.

	Paramters
	---------
	conf : ns_patches.config.get_conf() object
		Configuration info

	Returns
	-------
	trials : numpy array of ints, ( n_patches, n_trials )
		Each entry contains the image index for each patch and trial

	"""

	# generate a vector with the image indices repeated the number of repetitions
	img_seq = np.repeat( np.arange( conf.exp.n_img ),
	                     conf.exp.n_img_rep_per_run
	                   )

	# and shake them up
	np.random.shuffle( img_seq )

	# propogate the sequence for each of the patches
	trials = np.tile( img_seq, ( conf.stim.n_patches, 1 ) )

	# just check its how we think it is
	assert trials.shape == ( conf.stim.n_patches,
	                         conf.exp.n_img * conf.exp.n_img_rep_per_run
	                       )

	# calculate the sequence of incoherent patches per trial
	trial_incoh = sel_incoh( conf.exp.n_mod_patches,
	                         conf.exp.n_incoh_patches,
	                         conf.exp.n_trials
	                       )


	# we only need to consider the patches that will be modified; the others
	# already have their image set
	for i_incoh_patch in xrange( conf.exp.n_mod_patches ):

		# find out which patch, overall, this mod patch corresponds to
		i_patch = conf.exp.mod_patches[ i_incoh_patch ]

		# find the trials where this patch is designated as incoherent
		( _, incoh_trials ) = np.where( trial_incoh == i_incoh_patch )

		# when we do the shuffling of the incoherent events, we want to avoid
		# replacing an image with the same image, by chance. so we might need to
		# iterate a few times
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

	return trials


def test_sel2( n_patch, n_incoh, n_trials, n_img ):

	t = np.repeat( np.arange( n_img ), n_trials / n_img )

	t = np.tile( t, ( n_patch, 1 ) )

	success = False

	while not success:

		restart = False

		swapped = np.zeros( ( n_patch,  n_trials ) )
		swap_k = np.zeros( ( n_patch ) )

		n = 0

		while n < t.shape[ 1 ]:

			p = ( 1 - swap_k / 2.0 )

			p /= p.sum()

			try:
				s = np.random.choice( a = n_patch,
				                      size = n_incoh,
				                      replace = False,  # without replacement
				                      p = p
				                    )
			except ValueError:
				restart = True
				break

			for ( i_s, cs ) in enumerate( s ):

				p = ( 1 - np.double( swapped[ cs, : ] ) )
				p /= p.sum()

				try:
					c = np.random.choice( a = n_trials,
					                      size = 1,
					                      replace = False,  # without replacement
					                      p = p
					                    )
				except ValueError:
					restart = True
					break

				temp_1 = t[ cs, n ].copy()
				temp_2 = t[ cs, c ].copy()

				t[ cs, n ] = temp_2
				t[ cs, c ] = temp_1

				swapped[ cs, n ] = True
				swapped[ cs, c ] = True

				swap_k[ cs ] += 1

			n += 1


		if restart is not True:
			success = True


	return t






def test_sel( n_patch, n_incoh, n_trials, n_img ):

	success = False

	i = np.empty( ( n_patch, n_trials ) )

	while not success:

		img_k = np.zeros( ( n_img ) )
		patch_k = np.zeros( ( n_patch, n_img ) )

		restart = False

		for i_trial in xrange( n_trials ):

			for i_patch in xrange( n_patch ):

				p = 1 - ( img_k / 4.0 )
				p += 1 - ( np.sum( patch_k, axis= 0 ) / float( n_incoh ) )

				p[ p < 0 ] = 0.0

				p /= p.sum()

#				p[ np.isnan( p ) ] = 0.0

#				p = 1.0 - p

#				p /= p.sum()

				try:
					sel = np.random.choice( a = n_img,
					                        p = p
					                      )
				except ValueError:
					restart = True
					print p
					print ( i_trial, i_patch )
					break

				img_k[ sel ] += 1
				patch_k[ i_patch, sel ] += 1

				i[ i_patch, i_trial ] = sel

			if restart:
				break

		if not restart:
			success = True

	return i

def sel_incoh( n_patch, n_incoh, n_trials, img_seq ):
	"""Select the patches that will be incoherent in each trial.

	Parameters
	----------
	n_patch : integer
		Number of patches overall on a given trial
	n_incoh : integer
		Number of incoherent patches per trial
	n_trials : integer
		Total number of trials

	Returns
	-------
	incoh : numpy array of integers, ( incoh_patch, trial )
		Incoherent patch indices for each trial

	"""

	# allocate half the patches on each trial to be incoherent
	total_incoh = n_trials / 2

	# we may need to try a few times
	success = False

	while not success:

		restart = False

		# keep track of the number of incoherent events there has been in each
		# patch
		incoh_k = np.zeros( n_patch )

		incoh_img = np.zeros( ( n_patch, len( np.unique( img_seq ) ) ) )

		# initialise the container
		incoh = np.empty( ( n_incoh, n_trials ) )
		incoh.fill( np.NAN )

		for i_trial in xrange( n_trials ):

			trial_img = img_seq[ i_trial ]

			# want to make the probabilty of assigning a given patch to be incoherent
			# on this trial as proportional to a lack of assignment in previous
			# trials
			p = ( 1 - ( incoh_k / float( total_incoh ) ) )

			p += ( 1 - ( incoh_img[ :, trial_img ] / 2.0 ) )

			p /= p.sum()

			# `choice` throws a ValueError if there aren't enough non-zero
			# probability entries - meaning this sequence is a dud
			try:
				# choose a set of incoherent patches from a weighted sample of all
				# patches
				sel = np.random.choice( a = n_patch,
				                        size = n_incoh,
				                        replace = False,  # without replacement
				                        p = p
				                      )
			except ValueError:
				restart = True
				print i_trial
				break

			# assign the chosen patches
			incoh[ :, i_trial ] = sel

			# increment the counters for the chosen patches
			incoh_k[ sel ] += 1
			incoh_img[ sel, trial_img ] += 1

		# if we'd made it this far, and havent been asked to restart, then we have
		# a good set
		if not restart:
			success = True

	# ... but let's make sure
	# 1. each patch should have half the number of trials assigned as incoherent
	assert np.all( np.array( [ ( incoh == p ).sum() for p in xrange( n_patch ) ] )
	               == total_incoh
	             )

	# 2. pretty similar to the above
	assert np.all( np.histogram( incoh.flatten(), np.arange( n_patch + 1 ) )[ 0 ]
	               == total_incoh
	             )

	# 3. number of unique patches on each trial should equal the number of
	# incoherent patches
	assert np.all( np.array( [ len( np.unique( c ) ) for c in incoh.T ] )
	               == n_incoh
	             )

	return incoh
