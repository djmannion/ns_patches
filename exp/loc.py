import os.path
import csv

import numpy as np

import psychopy.visual, psychopy.misc, psychopy.event, psychopy.core

import ns_patches.config


def run( run_num ):

	conf = ns_patches.config.get_conf()

	win = psychopy.visual.Window( ( 1024, 768 ),
	                              monitor = "UMN_7T",
	                              fullscr = False,
	                              allowGUI = True
	                            )

	stim = get_stim( win, conf )

	run_clock = psychopy.core.Clock()

	quit_key = "q"
	trig_key = "t"

	run_clock.reset()

	if quit_key in keys:
		print "User aborted"
		win.close()
		return 1

	while run_clock.getTime() < conf.loc.full_run_len_s:

		# look up contrast

		# draw
		_ = [ patch.draw() for patch in stim ]

		win.flip()


	win.close()


def get_stim( win, conf ):

	stim = [ psychopy.visual.GratingStim( win = win,
	                                      tex = "sqrXsqr",
	                                      mask = "raisedCos",
	                                      pos = ( patch[ "px" ], patch[ "py" ] ),
	                                      units = "pix",
	                                      size = patch[ "diam" ],
	                                      sf = 5.0,
	                                      ori = 45.0,
	                                      phase = np.random.rand(),
	                                      contrast = 0.0
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




