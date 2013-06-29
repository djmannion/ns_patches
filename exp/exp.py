import os.path
import csv

import numpy as np

import psychopy.visual, psychopy.misc, psychopy.event, psychopy.core

import ns_patches.config, ns_patches.paths, ns_patches.prepare.prepare_images


def run( subj_id, run_num, show_perf = True ):

	# no matter what, can't show performance if its the first run
	if run_num == 1:
		show_perf = False

	conf = ns_patches.config.get_conf()

	paths = ns_patches.paths.get_exp_paths( conf )

	# this is an ( n_patches, n_trials ) array that contains the image info for
	# each patch and trial
	img_trials = ns_patches.config.gen_exp_patch_timing( conf )


	# this is an `n_run_seq` long vector of onset times, in s
	seq = gen_run_seq( conf )

	# work out how many of the 'pre' events contain stimuli
	n_stim_evts_in_pre = len( seq ) - conf.exp.n_trials

	# because we want to wraparound some images at the start
	img_trials = np.concatenate( ( img_trials[ :, -n_stim_evts_in_pre: ],
	                               img_trials
	                             ), axis = 1
	                           )

	assert img_trials.shape[ 1 ] == len( seq )

	win = psychopy.visual.Window( ( 1024, 768 ),
	                              monitor = "UMN_7T",
	                              fullscr = False,
	                              allowGUI = True
	                            )

	wait_str = "Loading stimuli..."

	wait_text = psychopy.visual.TextStim( win = win,
	                                      text = wait_str,
	                                      height = 22,
	                                      units = "pix",
	                                      bold = False,
	                                      pos = ( 0, -100 )
	                                    )

	wait_text.draw()
	win.flip()

	img = load_images( conf, paths )

	stim = ns_patches.prepare.prepare_images.set_stim( conf, win )

	stim_updated = [ False ] * len( seq )
	stim_drawn = [ False ] * len( seq )

	# prep for first trial
	stim = update_stim( stim, img, img_trials[ 0 ] )
	stim_updated[ 0 ] = True

	wait_str = "Press a button when ready for the run"
	wait_text.setText( wait_str )
	wait_text.draw()
	win.flip()

	keys = psychopy.event.waitKeys()

	wait_str = "Awaiting trigger..."
	wait_text.setText( wait_str )
	wait_text.draw()
	win.flip()

	run_clock = psychopy.core.Clock()

	quit_key = "q"
	trig_key = "t"

	keys = psychopy.event.waitKeys( keyList = [ quit_key, trig_key ] )

	run_clock.reset()

	if quit_key in keys:
		print "User aborted"
		win.close()
		return 1

	run_time = run_clock.getTime()

	while run_time < conf.exp.run_len_s:

		try:
			i_seq = np.where( run_time >= seq )[ 0 ][ -1 ]
			delta_t = run_time - seq[ i_seq ]
		except IndexError:
			delta_t = np.Inf

		if not stim_drawn[ i_seq ] and delta_t < conf.exp.img_on_s:

#			print delta_t

			map( lambda x : x.draw(), stim )
			win.flip()
			stim_drawn[ i_seq ] = True

			if i_seq < len( seq - 1 ):
				s = run_clock.getTime()
				stim = update_stim( stim, img, img_trials[ :, i_seq + 1 ] )
				print run_clock.getTime() - s

		if ( conf.exp.img_on_s < delta_t < conf.exp.bin_len_s ):

			win.flip()

			# look for a keypress
			keys = psychopy.event.getKeys()

			for key in keys:

				if key == quit_key:
					print "User abort"
					win.close()
					return 1

		run_time = run_clock.getTime()

	win.close()


def update_stim( stim, img, img_seq ):

	for i_stim in xrange( len( stim ) ):
		stim[ i_stim ].update_image( img[ img_seq[ i_stim ], ... ] )

	return stim



def load_images( conf, paths ):

	img_db_info = ns_patches.prepare.prepare_images.read_img_db_info( paths.img_db_info )

	i_img = np.where( img_db_info[ "status" ] == "Y" )[ 0 ]

	img_db = img_db_info[ i_img ]

#	assert len( img_db ) == conf.exp.n_img
	img_db = img_db[ :conf.exp.n_img ]

	img = np.empty( ( conf.exp.n_img,
	                  conf.stim.img_diam_pix,
	                  conf.stim.img_diam_pix
	                )
	              )
	img.fill( np.NAN )

	for ( i_curr, curr_img ) in enumerate( img_db ):

		img_path = os.path.join( paths.img_db,
		                         curr_img[ "album" ],
		                         curr_img[ "image" ]
		                       )

		img[ i_curr, ... ] = ns_patches.stimulus.load_img( img_path )

	assert np.sum( np.isnan( img ) ) == 0

	return img


def gen_run_seq( conf ):

	seq = np.ones( ( conf.exp.n_seq ) )
	seq[ :conf.exp.n_null ] = 0

	np.random.shuffle( seq )

	seq = np.concatenate( ( seq[ -conf.exp.n_pre: ], seq ) )

	assert len( seq ) == conf.exp.n_run_seq

	seq_s = np.where( seq != 0 )[ 0 ] * conf.exp.bin_len_s

	return seq_s


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



def set_stim( conf, stim, timing, run_time, ph_offs ):

	for ( patch, patch_t, ph_off ) in zip( stim, timing, ph_offs ):

		t_diff = run_time - np.array( patch_t )

		t_rect_diff = t_diff[ t_diff >= 0 ]

		if len( t_rect_diff ) == 0:
			contrast = 0.0
			continue

		t_off = np.min( t_rect_diff )

		if t_off <= conf.loc.dur_s:
			contrast = 1.0
		else:
			contrast = 0.0

		if ( np.mod( t_off + ph_off, conf.loc.reversal_interval_s * 2 ) <
		     conf.loc.reversal_interval_s
		   ):
			contrast *= -1.0

		patch.setContrast( contrast )

	return stim


def get_stim( win, conf ):

	stim = [ psychopy.visual.GratingStim( win = win,
	                                      tex = "sqrXsqr",
	                                      mask = conf.stim.mask_edge,
	                                      pos = ( patch[ "cx" ], patch[ "cy" ] ),
	                                      units = "pix",
	                                      size = patch[ "diam" ],
	                                      sf = 5.0,
	                                      ori = 45.0,
	                                      phase = np.random.rand(),
	                                      contrast = 0.0,
	                                      maskParams = { "fringeWidth" :
	                                                     conf.stim.mask_edge_frac
	                                                   }
	                                    )
	        for patch in conf.stim.patches
	      ]

	return stim


def load_timing( conf, exp_paths, run_num ):

	timing = []

	for i_patch in xrange( conf.stim.n_patches ):

		with open( os.path.join( exp_paths.timing_dir,
		                         ( "ns_patches-loc_timing_patch_" +
		                           "{n:02d}.txt".format( n = i_patch )
		                         )
		                       ), "r"
		         ) as timing_file:

			timing_csv = csv.reader( timing_file, delimiter = " " )

			p_timing = [ run_timing for run_timing in timing_csv ]

			timing.append( map( int, p_timing[ run_num - 1 ] ) )

	return timing


def init_task( conf ):
	"""Initialises the task timing.

	Returns
	-------
	task_lut : numpy array, shape of ( evt x info )
		Task lookup table, where dim two is ( time_s, digit, polarity, target )
	targets : numpy array, shape of ( target, info )
		Target information, stored as ( digit, polarity )

	"""

	n_task_per_run = int( conf.loc.run_len_s * 3.0 )

	task_set = np.arange( 10 )
	np.random.shuffle( task_set )

	targets = np.array( [ [ task_set[ i ], [ -1, +1 ][ i ] ]
	                      for i in xrange( 2 )
	                    ]
	                  )

	# second dim is (time, digit, polarity, target or not)
	task_lut = np.empty( ( n_task_per_run, 4 ) )

	for i_evt in xrange( n_task_per_run ):

		time_s = i_evt * ( 1.0 / 3.0 )

		curr_task_set = task_set.copy()
		curr_task_set = curr_task_set[ curr_task_set != task_lut[ i_evt - 1, 1 ] ]

		digit = curr_task_set[ np.random.randint( len( curr_task_set ) ) ]

		polarity = [ -1, +1 ][ np.random.randint( 2 ) ]

		if np.any( np.logical_and( targets[ :, 0 ] == digit,
		                           targets[ :, 1 ] == polarity
		                         )
		         ):
			target = 1
		else:
			target = 0

		task_lut[ i_evt, : ] = [ time_s, digit, polarity, target ]

	return ( task_lut, targets )

