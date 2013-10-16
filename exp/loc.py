import os.path
import csv

import numpy as np

import psychopy.visual, psychopy.misc, psychopy.event, psychopy.core

import ns_patches.config, ns_patches.paths, ns_patches.stimulus


def run( run_num ):

    conf = ns_patches.config.get_conf()

    paths = ns_patches.paths.get_exp_paths( conf )

    timing = load_timing( conf, paths, run_num )

    win = psychopy.visual.Window( ( 1024, 768 ),
                                  monitor = "UMN_7T_colour",
                                  fullscr = False,
                                  allowGUI = True
                                )

    stim = get_stim( win, conf )

    ph_offs = np.random.rand( len( stim ) ) * 0.15

    fix_stim = ns_patches.stimulus.get_fixation( win, conf )

    ( task, targs ) = init_task( conf )

    fix_text = psychopy.visual.TextStim( win = win,
                                         text = "",
                                         height = 26,
                                         units = "pix",
                                         bold = False
                                       )

    targ_pos = ( -62, 62 )

    targ_txt = [ psychopy.visual.TextStim( win = win,
                                           text = "{s:d}".format( s = targs[ i, 0 ] ),
                                           height = 26,
                                           units = "pix",
                                           pos = ( targ_pos[ i ], 0 ),
                                           color = np.repeat( targs[ i, 1 ], 3 )
                                         )
                 for i in xrange( 2 )
               ]


    wait_str = "Press a button when ready for the next run"

    wait_text = psychopy.visual.TextStim( win = win,
                                          text = wait_str,
                                          height = 22,
                                          units = "pix",
                                          bold = False,
                                          pos = ( 0, -100 )
                                        )


    run_clock = psychopy.core.Clock()

    quit_key = "q"
    trig_key = "t"

    _ = [ fixation.draw() for fixation in fix_stim ]
    _ = [ t_text.draw() for t_text in targ_txt ]

    wait_text.draw()
    win.flip()

    win.getMovieFrame()
    win.saveMovieFrames( "caps/loc.png" )

    keys = psychopy.event.waitKeys()

    _ = [ fixation.draw() for fixation in fix_stim ]
    _ = [ t_text.draw() for t_text in targ_txt ]

    win.flip()

    keys = psychopy.event.waitKeys( keyList = [ quit_key, trig_key ] )

    run_clock.reset()

    if quit_key in keys:
        print "User aborted"
        win.close()
        return 1

    run_time = run_clock.getTime()

    while run_time < conf.loc.run_len_s:

        stim = set_stim( conf, stim, timing, run_time + ( 1.0 / 60.0 ), ph_offs )

        # draw
        _ = [ fixation.draw() for fixation in fix_stim[ :-1 ] ]
        _ = [ patch.draw() for patch in stim ]

        i_task_evt = np.where( run_time > task[ :, 0 ] )[ 0 ][ -1 ]

        fix_text.setText( str( int( task[ i_task_evt, 1 ] ) ) )
        fix_text.setColor( np.repeat( task[ i_task_evt, 2 ], 3 ) )

        fix_text.draw()

        fix_stim[ -1 ].draw()

        win.flip()
        run_time = run_clock.getTime()

        keys = psychopy.event.getKeys()

        for key in keys:

            if key == quit_key:
                print "User abort"
                win.close()
                return 1

    win.close()


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

