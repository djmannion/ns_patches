import os.path

import numpy as np
import scipy.stats

import psychopy.visual, psychopy.misc, psychopy.event, psychopy.core

import ns_patches.config, ns_patches.paths, ns_patches.prepare.prepare_images
import ns_patches.exp.loc

def run( subj_id, run_num, show_perf = True, mon_name = "UMN_7T_colour" ):

	# no matter what, can't show performance if its the first run
	if run_num == 1:
		show_perf = False

	conf = ns_patches.config.get_conf()

	paths = ns_patches.paths.get_exp_paths( conf )

	# get the paths for the relevant output, making sure they don't exist already
	log_path = get_log_path( subj_id, run_num, paths )

	# if we want to show the last run performance, need to load its output
	if show_perf:

		prev_log_path = get_log_path( subj_id,
		                              run_num - 1,
		                              paths,
		                              err_if_exist = False
		                            )

		prev_perf = eval_resp( np.load( prev_log_path )[ "resp_status" ] )

		perf_str = "Prev run: " + prev_perf[ 1 ]

	else:
		perf_str = ""


	# this is an ( n_patches, n_trials ) array that contains the image info for
	# each patch and trial
	img_trials = gen_patch_timing( conf )

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

	resp_status = np.empty( len( seq ) )
	resp_status.fill( np.NAN )

	win = psychopy.visual.Window( ( 1024, 768 ),
	                              monitor = mon_name,
	                              fullscr = True,
	                              allowGUI = False
	                            )

	perf_text = psychopy.visual.TextStim( win = win,
	                                      text = perf_str,
	                                      height = 22,
	                                      units = "pix",
	                                      bold = False,
	                                      pos = ( 0, -150 )
	                                    )

	fix_stim = ns_patches.stimulus.get_fixation( win, conf )

	wait_str = "Loading stimuli..."

	wait_text = psychopy.visual.TextStim( win = win,
	                                      text = wait_str,
	                                      height = 22,
	                                      units = "pix",
	                                      bold = False,
	                                      pos = ( 0, -100 )
	                                    )

	map( lambda x : x.draw(), fix_stim )
	wait_text.draw()
	perf_text.draw()
	win.flip()

	img_db = load_img_info( conf, paths )
	img = load_images( conf, paths, img_db )
	masks = gen_masks( conf )
	stim = gen_stim( conf, win, masks )

	stim_updated = [ False ] * len( seq )
	stim_drawn = [ False ] * len( seq )

	# prep for first trial
	stim = update_stim( conf, stim, img, masks, img_trials[ :, 0 ] )
	stim_updated[ 0 ] = True

	map( lambda x : x.draw(), fix_stim )
	wait_str = "Press a button when ready for the run"
	wait_text.setText( wait_str )
	wait_text.draw()
	perf_text.draw()
	win.flip()

	keys = psychopy.event.waitKeys()

	map( lambda x : x.draw(), fix_stim )
	wait_str = "Awaiting trigger..."
	wait_text.setText( wait_str )
	wait_text.draw()
	perf_text.draw()
	win.flip()

	wait_text.setPos( [ 0, -350 ] )

	run_clock = psychopy.core.Clock()

	quit_key = "q"
	trig_key = "t"
	resp_keys = [ "b", "r" ]

	keys = psychopy.event.waitKeys( keyList = [ quit_key, trig_key ] )

	run_clock.reset()

	if quit_key in keys:
		print "User aborted"
		win.close()
		return 1

	wait_text.setText( "" )
	map( lambda x : x.draw(), fix_stim )
	wait_text.draw()
	win.flip()

	run_time = run_clock.getTime()

	while run_time < conf.exp.run_len_s:

		# handle a situation where the first n events might be nulls
		if run_time < seq[ 0 ]:
			delta_t = np.Inf
		else:
			# find the sequence index for this timepoint
			i_seq = np.where( run_time >= seq )[ 0 ][ -1 ]
			delta_t = run_time - seq[ i_seq ]

		# if it is the stimulus display period, but a stimulus hasn't actually been
		# displayed
		if delta_t < conf.exp.img_on_s and not stim_drawn[ i_seq ]:

			map( lambda x : x.draw(), fix_stim )
			stim.draw()
			wait_text.draw()
			win.flip()
			stim_drawn[ i_seq ] = True

			# we don't want any previous keypresses hanging around
			psychopy.event.clearEvents()

			# find the index of the coherent image currently being shown
			i_img = int( scipy.stats.mode( img_trials[ :, i_seq ] )[ 0 ] )
			img_album = img_db[ i_img ][ "album" ]

			win.getMovieFrame()
			win.saveMovieFrames( "caps/capn{t:02d}.png".format( t = i_seq ) )

			# take this time to update our next stimulus
			if i_seq < ( len( seq ) - 1 ):
				stim = update_stim( conf, stim, img, masks, img_trials[ :, i_seq + 1 ] )
				stim_updated[ i_seq + 1 ] = True

		# we've now entered the stimulus off period, but still within the limits of
		# this bin. here, we want to gather any responses from the subject
		if ( conf.exp.img_on_s < delta_t < conf.exp.bin_len_s ):

			# look for a keypress. this also returns any keys that were pressed
			# during the image display period
			keys = psychopy.event.getKeys( keyList = resp_keys + [ quit_key ] )

			for key in keys:

				if key == quit_key:
					print "User abort"
					win.close()
					return 1

				# not quit, and since we were only listening for it and the response
				# keys, it must be a response
				else:

					corr = np.NAN

					# work through the possibilities
					if img_album == "Flowers":
						if key == "r":
							corr = 1
						else:
							corr = 0
					else:
						if key == "b":
							corr = 1
						else:
							corr = 0

					# update the response status with the response correctness
					resp_status[ i_seq ] = corr

			if i_seq > 0:
				resp_stat = eval_resp( resp_status[ :( i_seq - 1 ) ] )
				wait_text.setText( resp_stat[ 1 ] )

			map( lambda x : x.draw(), fix_stim )

			wait_text.draw()

			win.flip()

		run_time = run_clock.getTime()

	win.close()

	assert np.all( np.array( stim_drawn ) == True )
	assert np.all( np.array( stim_updated ) == True )

	np.savez( file = log_path,
	          img_trials = img_trials,
	          seq = seq,
	          resp_status = resp_status
	        )


def eval_resp( resp_status ):

	n_corr = np.nansum( resp_status )
	n_tot = np.sum( np.logical_not( np.isnan( resp_status ) ) )
	n_miss = np.sum( np.isnan( resp_status ) )

	if np.isnan( n_corr ):
		n_corr = 0

	resp_text = "{c:.0f} / {t:.0f}, {n:.0f} missed".format( c = n_corr,
	                                                        t = n_tot,
	                                                        n = n_miss
	                                                      )

	return ( ( n_corr, n_tot, n_miss ), resp_text )


def get_log_path( subj_id, run_num, paths, err_if_exist = True ):

	log_file = "{s:s}_ns_patches-run_{r:d}_log.npz".format( s = subj_id,
	                                                        r = run_num,
	                                                      )

	log_path = os.path.join( paths.log_dir, log_file )

	if os.path.exists( log_path ) and err_if_exist:
		raise IOError( "Output path " + log_path + " already exists" )

	return log_path



def update_stim( conf, stim, img, masks, img_seq ):

	new_img = np.zeros( [ conf.stim.img_diam_pix ] * 2 + [ 3 ] )

	for i_patch in xrange( len( img_seq ) ):

		i_img = img_seq[ i_patch ]

		for i_chan in xrange( 3 ):
			new_img[ ..., i_chan ] += ( img[ i_img, :, :, i_chan ] *
			                            ( masks[ i_patch, ... ] > 0 )
			                          )

	stim.setImage( new_img )

	return stim


def gen_patch_timing( conf ):

	# this is 80 x 26
	mod_design = ns_patches.config.gen_exp_patch_timing( conf )

	i_mod_patches = conf.exp.mod_patches
	i_unmod_patches = np.setdiff1d( range( conf.stim.n_patches ),
	                                i_mod_patches
	                              )

	design = np.empty( ( conf.stim.n_patches, conf.exp.n_trials ) )
	design.fill( np.NAN )

	for i_trial in xrange( conf.exp.n_trials ):

		mod_patches = mod_design[ i_trial, : ]

		trial_coh = scipy.stats.mode( mod_patches )[ 0 ]

		design[ i_unmod_patches, i_trial ] = trial_coh
		design[ i_mod_patches, i_trial ] = mod_design[ i_trial, : ]

	return design


def gen_stim( conf, win, masks ):

	img_size = [ conf.stim.img_diam_pix ] * 2

	img_mask = np.sum( masks, axis = 0 ) * 2.0 - 1.0

	stim = psychopy.visual.ImageStim( win = win,
	                                  image = np.ones( img_size ),
	                                  size = img_size,
	                                  units = "pix",
	                                  mask = img_mask
	                                )

	return stim


def gen_masks( conf ):

	masks = np.empty( ( conf.stim.n_patches,
	                    conf.stim.img_diam_pix,
	                    conf.stim.img_diam_pix
	                  )
	                )
	masks.fill( np.NAN )

	mask_f = ns_patches.stimulus.make_mask

	for ( i_mask, patch ) in enumerate( conf.stim.patches ):

		masks[ i_mask, ... ] = mask_f( patch = patch,
		                               img_size = masks.shape[ 1 ],
		                               mask_range = [ 0, 1 ]
		                             )

	return masks


def load_img_info( conf, paths ):

	img_db_info = ns_patches.prepare.prepare_images.read_img_db_info( paths.img_db_info )

	i_img = np.where( img_db_info[ "status" ] == "Y" )[ 0 ]

	img_db = img_db_info[ i_img ]

	assert len( img_db ) == conf.exp.n_img

	return img_db


def load_images( conf, paths, img_db ):

	img = np.empty( ( conf.exp.n_img,
	                  conf.stim.img_diam_pix,
	                  conf.stim.img_diam_pix,
	                  3
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


def check_run_seq( conf ):

	seq = gen_run_seq( conf )

	cmd = [ "3dDeconvolve",
	        "-nodata",
	        "{n:0f}".format( n = conf.exp.n_vol ),
	        "{n:0f}".format( n = conf.acq.tr_s ),
	        "-polort", "A",
	        "-CENSORTR", "0-{n:0f}".format( n = conf.exp.n_censor_vols - 1 ),
	        "-num_stimts", "1",
	        "-local_times",
	        "-stim_times_IM",
	        "1",
	        "'1D: " + " ".join( map( str, seq ) ) + "'",
	        "'SPMG1(1)'"
	      ]

	print " ".join( cmd )

