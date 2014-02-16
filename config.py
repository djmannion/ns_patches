"Configuration file for the natural scenes patches fMRI study"

import os.path
import csv

import numpy as np
import scipy.stats

import psychopy.misc

import ns_patches.paths


class ConfigContainer( object ):
    pass


def get_conf( subj_id = None, subj_types = "all", excl_subj_ids = None ):

    if excl_subj_ids is None:
        excl_subj_ids = []

    conf = ConfigContainer()

    conf.stim = _get_stim_conf()
    conf.acq = _get_acq_conf()
    conf.exp = _get_exp_conf( conf )
    conf.loc = _get_loc_conf( conf )
    conf.ana = _get_ana_conf()
    conf.all_subj = _get_subj_conf( None, subj_types, excl_subj_ids )

    if subj_id is not None:
        conf.subj = _get_subj_conf( subj_id )

    return conf


def _get_exp_conf( conf ):

    exp_conf = ConfigContainer()

    exp_conf.id = "ns_patches"

    exp_conf.mod_patches = range( 28 )
    # remove patches on the vertical meridian
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

    ana_conf.exclude_subj_ids = [ "s1023" ]  # poor localisers

    # where( any( k < 10, axis = 0 ) ), after getting rid of ^
    ana_conf.exclude_patch_ids = [ 8, 18, 20, 21, 26 ]

    return ana_conf


def _get_acq_conf():

    acq_conf = ConfigContainer()

    acq_conf.monitor_name = "UMN_7T_colour"
    acq_conf.test_monitor_name = "N13_CRT"

    acq_conf.tr_s = 2.0

    # how to reshape the data to be in +RAS convention
    # subsequent commands are relative to the data AFTER this operation
    # see docs for how to determine these
    acq_conf.ras = ( "-x", "-z", "-y" )

    # phase encode direction, according to the data's internal axes
    acq_conf.ph_enc_dir = "x"

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

    # t-threshold for analysis
    loc_conf.t_p = 0.001

    # degrees of freedom
    loc_conf.dof = 844

    loc_conf.area_thr = 5

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

    # ~ 1.8, 3.5, 6.1, 9.75 deg
    ring_ecc = [ 58, 113, 196, 313 ]

    # M0 = 29.2/3.67; M = 29.2/(e+3.67); D = M0/M*15*2*1.08
    ring_diam = [ 45, 59, 80, 113 ]
    ring_diam = map( lambda x : x * ( 1 / 1.08 ), ring_diam )

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


def _get_subj_conf( subj_id = None, subj_types = "all", excl_subj_ids = None ):

    if excl_subj_ids is None:
        excl_subj_ids = []

    s1000 = ConfigContainer()

    s1000.subj_id = "s1000"
    s1000.fs_subj_id = "s1000"
    s1000.acq_date = "20130813"
    s1000.comments = ""
    s1000.n_runs = 8
    s1000.mot_base = 5
    s1000.vol_base = 108
    s1000.is_loc = False
    s1000.mask_SI = 90

    s1000.extra_al_params = [ "-parang", "1", "-12", "-2",
                              "-parang", "2", "14", "24",
                              "-parang", "3", "20", "37",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1000.node_k = { "lh" : 130318,
                     "rh" : 131151
                   }

    s1000_loc = ConfigContainer()

    s1000_loc.subj_id = "s1000_loc"
    s1000_loc.fs_subj_id = "s1000"
    s1000_loc.acq_date = "20130812"
    s1000_loc.comments = ""
    s1000_loc.n_runs = 6
    s1000_loc.mot_base = 4
    s1000_loc.vol_base = 83
    s1000_loc.is_loc = True
    s1000_loc.mask_SI = 85

    s1000_loc.extra_al_params = [ "-parang", "1", "-4", "6",
                                  "-parang", "2", "14", "24",
                                  "-parang", "3", "15", "32",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1000_loc.node_k = { "lh" : 130318,
                         "rh" : 131151
                       }

    s1021_loc = ConfigContainer()

    s1021_loc.subj_id = "s1021_loc"
    s1021_loc.fs_subj_id = "s1021"
    s1021_loc.acq_date = "20130812"
    s1021_loc.comments = ""
    s1021_loc.n_runs = 6
    s1021_loc.mot_base = 4
    s1021_loc.vol_base = 83
    s1021_loc.is_loc = True
    s1021_loc.mask_SI = 75

    s1021_loc.extra_al_params = [ "-parang", "1", "-5", "5",
                                  "-parang", "2", "11", "21",
                                  "-parang", "3", "46", "62",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1021_loc.node_k = { "lh" : 140847,
                         "rh" : 141381
                       }

    s1021 = ConfigContainer()

    s1021.subj_id = "s1021"
    s1021.fs_subj_id = "s1021"
    s1021.acq_date = "20130813"
    s1021.comments = ""
    s1021.n_runs = 8
    s1021.mot_base = 5
    s1021.vol_base = 108
    s1021.is_loc = False
    s1021.mask_SI = 75

    s1021.extra_al_params = [ "-parang", "1", "-8", "2",
                              "-parang", "2", "13", "23",
                              "-parang", "3", "46", "56",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1021.node_k = { "lh" : 140847,
                     "rh" : 141381
                   }

    s1008 = ConfigContainer()

    s1008.subj_id = "s1008"
    s1008.fs_subj_id = "s1008"
    s1008.acq_date = "20130823"
    s1008.comments = ""
    s1008.n_runs = 8
    s1008.mot_base = 5
    s1008.vol_base = 108
    s1008.is_loc = False
    s1008.mask_SI = 100

    s1008.extra_al_params = [ "-parang", "1", "-5", "5",
                              "-parang", "2", "5", "15",
                              "-parang", "3", "6", "17",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1008.node_k = { "lh" : 140427,
                     "rh" : 141898
                   }


    s1008_loc = ConfigContainer()

    s1008_loc.subj_id = "s1008_loc"
    s1008_loc.fs_subj_id = "s1008"
    s1008_loc.acq_date = "20130822"
    s1008_loc.comments = ""
    s1008_loc.n_runs = 6
    s1008_loc.mot_base = 4
    s1008_loc.vol_base = 83
    s1008_loc.is_loc = True
    s1008_loc.mask_SI = 100

    s1008_loc.extra_al_params = [ "-parang", "1", "-5", "7",
                                  "-parang", "2", "9", "19",
                                  "-parang", "3", "6", "16",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1008_loc.node_k = { "lh" : 140427,
                         "rh" : 141898
                       }

    s1023 = ConfigContainer()

    s1023.subj_id = "s1023"
    s1023.fs_subj_id = "s1023"
    s1023.acq_date = "20130823"
    s1023.comments = ""
    s1023.n_runs = 10
    s1023.mot_base = 5
    s1023.vol_base = 108
    s1023.is_loc = False
    s1023.mask_SI = 95

    s1023.extra_al_params = [ "-parang", "1", "-4", "6",
                              "-parang", "2", "16", "26",
                              "-parang", "3", "45", "58",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1023.node_k = { "lh" : 132517,
                     "rh" : 130499
                   }

    s1023_loc = ConfigContainer()

    s1023_loc.subj_id = "s1023_loc"
    s1023_loc.fs_subj_id = "s1023"
    s1023_loc.acq_date = "20130822"
    s1023_loc.comments = ""
    s1023_loc.n_runs = 6
    s1023_loc.mot_base = 4
    s1023_loc.vol_base = 83
    s1023_loc.is_loc = True
    s1023_loc.mask_SI = 95

    s1023_loc.extra_al_params = [ "-parang", "1", "-5", "5",
                                  "-parang", "2", "24", "34",
                                  "-parang", "3", "30", "40",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1023_loc.node_k = { "lh" : 132517,
                         "rh" : 130499
                       }

    s1011 = ConfigContainer()

    s1011.subj_id = "s1011"
    s1011.fs_subj_id = "s1011"
    s1011.acq_date = "20130828"
    s1011.comments = ""
    s1011.n_runs = 8
    s1011.mot_base = 5
    s1011.vol_base = 108
    s1011.is_loc = False
    s1011.mask_SI = 100

    s1011.extra_al_params = [ "-parang", "1", "-5", "9",
                              "-parang", "2", "4", "14",
                              "-parang", "3", "10", "20",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1011.node_k = { "lh" : 128434,
                     "rh" : 128461
                   }

    s1011_loc = ConfigContainer()

    s1011_loc.subj_id = "s1011_loc"
    s1011_loc.fs_subj_id = "s1011"
    s1011_loc.acq_date = "20130827"
    s1011_loc.comments = ""
    s1011_loc.n_runs = 6
    s1011_loc.mot_base = 4
    s1011_loc.vol_base = 83
    s1011_loc.is_loc = True
    s1011_loc.mask_SI = 100

    s1011_loc.extra_al_params = [ "-parang", "1", "0", "13",
                                  "-parang", "2", "5", "15",
                                  "-parang", "3", "11", "17",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1011_loc.node_k = { "lh" : 128434,
                         "rh" : 128461
                       }

    s1048 = ConfigContainer()

    s1048.subj_id = "s1048"
    s1048.fs_subj_id = "s1048"
    s1048.acq_date = "20130829"
    s1048.comments = ""
    s1048.n_runs = 8
    s1048.mot_base = 5
    s1048.vol_base = 108
    s1048.is_loc = False
    s1048.mask_SI = 100

    s1048.extra_al_params = [ "-parang", "1", "-9", "5",
                              "-parang", "2", "5", "15",
                              "-parang", "3", "20", "34",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1048.node_k = { "lh" : 141613,
                     "rh" : 141576
                   }

    s1048_loc = ConfigContainer()

    s1048_loc.subj_id = "s1048_loc"
    s1048_loc.fs_subj_id = "s1048"
    s1048_loc.acq_date = "20130828"
    s1048_loc.comments = ""
    s1048_loc.n_runs = 6
    s1048_loc.mot_base = 4
    s1048_loc.vol_base = 83
    s1048_loc.is_loc = True
    s1048_loc.mask_SI = 100

    s1048_loc.extra_al_params = [ "-parang", "1", "-9", "5",
                                  "-parang", "2", "5", "18",
                                  "-parang", "3", "29", "41",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1048_loc.node_k = { "lh" : 141613,
                         "rh" : 141576
                       }

    s1046 = ConfigContainer()

    s1046.subj_id = "s1046"
    s1046.fs_subj_id = "s1046"
    s1046.acq_date = "20130911"
    s1046.comments = ""
    s1046.n_runs = 8
    s1046.mot_base = 4
    s1046.vol_base = 108
    s1046.is_loc = False
    s1046.mask_SI = 100

    s1046.extra_al_params = [ "-parini", "1", "-3",
                              "-parini", "2", "38",
                              "-parini", "3", "35",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1046.node_k = { "lh" : 127459,
                     "rh" : 125182
                   }


    s1046_loc = ConfigContainer()

    s1046_loc.subj_id = "s1046_loc"
    s1046_loc.fs_subj_id = "s1046"
    s1046_loc.acq_date = "20131011"  #"20130906"
    s1046_loc.comments = ""
    s1046_loc.n_runs = 6
    s1046_loc.mot_base = 4
    s1046_loc.vol_base = 83
    s1046_loc.is_loc = True
    s1046_loc.mask_SI = 100

    s1046_loc.extra_al_params = [ "-parini", "1", "4",
                                  "-parini", "2", "35",
                                  "-parini", "3", "6",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1046_loc.node_k = { "lh" : 127459,
                         "rh" : 125182
                       }

    s1033 = ConfigContainer()

    s1033.subj_id = "s1033"
    s1033.fs_subj_id = "s1033"
    s1033.acq_date = "20130916"
    s1033.comments = ""
    s1033.n_runs = 8
    s1033.mot_base = 5
    s1033.vol_base = 108
    s1033.is_loc = False
    s1033.mask_SI = 100

    s1033.extra_al_params = [ "-parini", "1", "5",
                              "-parini", "2", "29",
                              "-parini", "3", "40",
                              "-parini", "6", "-5",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1033.node_k = { "lh" : 144028,
                     "rh" : 145055
                   }



    s1033_loc = ConfigContainer()

    s1033_loc.subj_id = "s1033_loc"
    s1033_loc.fs_subj_id = "s1033"
    s1033_loc.acq_date = "20130912"
    s1033_loc.comments = ""
    s1033_loc.n_runs = 6
    s1033_loc.mot_base = 4
    s1033_loc.vol_base = 83
    s1033_loc.is_loc = True
    s1033_loc.mask_SI = 100

    s1033_loc.extra_al_params = [ "-parini", "1", "13",# "18",
                                  "-parini", "2", "26", #"31",
                                  "-parini", "3", "50",# "55",
                                  "-parini", "6", "-5",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1033_loc.node_k = { "lh" : 144028,
                         "rh" : 145055
                       }

    s1012 = ConfigContainer()

    s1012.subj_id = "s1012"
    s1012.fs_subj_id = "s1012"
    s1012.acq_date = "20130920"
    s1012.comments = ""
    s1012.n_runs = 8
    s1012.mot_base = 5
    s1012.vol_base = 108
    s1012.is_loc = False
    s1012.mask_SI = 100

    s1012.extra_al_params = [ "-parini", "1", "-2",# "5",
                              "-parini", "2", "35",# "43",
                              "-parini", "3", "40",# "27",
                              "-parini", "6", "-5",
                              "-maxrot", "10",
                              "-source_automask+2",
                              "-nocmass"
                            ]

    s1012.node_k = { "lh" : 155534,
                     "rh" : 159875
                   }

    s1012_loc = ConfigContainer()

    s1012_loc.subj_id = "s1012_loc"
    s1012_loc.fs_subj_id = "s1012"
    s1012_loc.acq_date = "20130917"
    s1012_loc.comments = ""
    s1012_loc.n_runs = 6
    s1012_loc.mot_base = 4
    s1012_loc.vol_base = 83
    s1012_loc.is_loc = True
    s1012_loc.mask_SI = 100

    s1012_loc.extra_al_params = [ "-parini", "1", "0",# "5",
                                  "-parini", "2", "38",# "43",
                                  "-parini", "3", "22",# "27",
                                  "-parini", "6", "-5",
                                  "-maxrot", "10",
                                  "-source_automask+2",
                                  "-nocmass"
                                ]

    s1012_loc.node_k = { "lh" : 155534,
                         "rh" : 159875
                       }


    subj = ConfigContainer()

    subj.subj = { "s1000" : s1000, "s1000_loc" : s1000_loc,
                  "s1021" : s1021, "s1021_loc" : s1021_loc,
                  "s1008" : s1008, "s1008_loc" : s1008_loc,
                  "s1023" : s1023, "s1023_loc" : s1023_loc,
                  "s1011" : s1011, "s1011_loc" : s1011_loc,
                  "s1048" : s1048, "s1048_loc" : s1048_loc,
                  "s1046" : s1046, "s1046_loc" : s1046_loc,
                  "s1033" : s1033, "s1033_loc" : s1033_loc,
                  "s1012" : s1012, "s1012_loc" : s1012_loc
                }

    if subj_types == "loc":
        subj.subj = { x : subj.subj[ x ] for x in subj.subj if "loc" in x }
    elif subj_types == "exp":
        subj.subj = { x : subj.subj[ x ] for x in subj.subj if not "loc" in x }

    for excl_subj_id in excl_subj_ids:
        del( subj.subj[ excl_subj_id ] )

    if subj_id is None:
        return subj
    else:
        return subj.subj[ subj_id ]

