import os.path

import numpy as np
import scipy.stats

import psychopy.visual, psychopy.misc, psychopy.event, psychopy.core

import ns_patches.config, ns_patches.paths, ns_patches.prepare.prepare_images
import ns_patches.exp.loc

def run( subj_id, run_num, show_perf = True, mon_name = "UMN_7T" ):

	# no matter what, can't show performance if its the first run
	if run_num == 1:
		show_perf = False

	conf = ns_patches.config.get_conf()

	paths = ns_patches.paths.get_exp_paths( conf )

	# get the paths for the relevant output, making sure they don't exist already
	log_path = get_log_path( subj_id, run_num, paths )

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

	resp_status = np.empty( len( seq ) )
	resp_status.fill( np.NAN )

	win = psychopy.visual.Window( ( 1024, 768 ),
	                              monitor = mon_name,
	                              fullscr = True,
	                              allowGUI = False
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
	win.flip()

	keys = psychopy.event.waitKeys()

	map( lambda x : x.draw(), fix_stim )
	wait_str = "Awaiting trigger..."
	wait_text.setText( wait_str )
	wait_text.draw()
	win.flip()

	wait_text.setPos( [ 0, -300 ] )

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

	run_time = run_clock.getTime()

	while run_time < conf.exp.run_len_s:

		if run_time < seq[ 0 ]:
			delta_t = np.Inf
		else:
			i_seq = np.where( run_time >= seq )[ 0 ][ -1 ]
			delta_t = run_time - seq[ i_seq ]

		if delta_t < conf.exp.img_on_s and not stim_drawn[ i_seq ]:

			map( lambda x : x.draw(), fix_stim )
			stim.draw()
			win.flip()
			stim_drawn[ i_seq ] = True


			psychopy.event.clearEvents()

			# find the index of the coherent image currently being shown
			i_img = int( scipy.stats.mode( img_trials[ :, i_seq ] )[ 0 ] )
			img_album = img_db[ i_img ][ "album" ]

#			win.getMovieFrame()
#			win.saveMovieFrames( "caps/cap{t:02d}.png".format( t = i_seq ) )

			if i_seq < ( len( seq ) - 1 ):
				stim = update_stim( conf, stim, img, masks, img_trials[ :, i_seq + 1 ] )
				stim_updated[ i_seq + 1 ] = True

		if ( conf.exp.img_on_s < delta_t < conf.exp.bin_len_s ):

			# look for a keypress
			keys = psychopy.event.getKeys( keyList = resp_keys + [ quit_key ] )

			for key in keys:

				if key == quit_key:
					print "User abort"
					win.close()
					return 1

				else:

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

					resp_status[ i_seq ] = corr

			wait_text.setText( "test" + str( resp_status[ i_seq ] ) )

			map( lambda x : x.draw(), fix_stim )

			wait_text.draw()

			win.flip()

		run_time = run_clock.getTime()

	win.close()

	assert np.all( np.array( stim_drawn ) == True )
	assert np.all( np.array( stim_updated ) == True )


def get_log_path( subj_id, run_num, paths, err_if_exist = True ):

	log_file = "{s:s}_ns_patches-run_{r:d}_log.npz".format( s = subj_id,
	                                                        r = run_num,
	                                                      )

	log_path = os.path.join( paths.log_dir, log_file )

	if os.path.exists( log_path ):
		raise IOError( "Output path " + log_path + " already exists" )

	return log_path



def update_stim( conf, stim, img, masks, img_seq ):

	new_img = np.zeros( [ conf.stim.img_diam_pix ] * 2 )

	for i_patch in xrange( len( img_seq ) ):

		i_img = img_seq[ i_patch ]

		new_img += img[ i_img, ... ] * ( masks[ i_patch, ... ] > 0 )

	stim.setImage( new_img )

	return stim


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

#	assert len( img_db ) == conf.exp.n_img
	img_db = img_db[ :conf.exp.n_img ]

	return img_db


def load_images( conf, paths, img_db ):

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
