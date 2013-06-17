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
