
import numpy as np
import scipy.stats
import scikits.bootstrap

import ns_patches.config, ns_patches.paths


def loc_stats( conf, paths ):

    all_conf = ns_patches.config.get_conf( subj_id = None, subj_types = "loc" )

    subj_ids = all_conf.all_subj.subj.keys()
    subj_ids.sort()

    patch_k = np.empty( ( len( subj_ids ), conf.exp.n_mod_patches ) )
    patch_k.fill( np.NAN )

    for ( i_subj, subj_id ) in enumerate( subj_ids ):

        subj_conf = ns_patches.config.get_conf( subj_id )
        subj_paths = ns_patches.paths.get_subj_paths( subj_conf )

        subj_patch_k = np.loadtxt( subj_paths.loc.patch_id_count.full( ".txt" ) )

        patch_k[ i_subj, : ] = subj_patch_k

    assert np.sum( np.isnan( patch_k ) ) == 0

    np.savetxt( paths.patch_k.full( ".txt" ), patch_k, fmt = "%d" )


def coh_summ( conf, paths ):

    all_conf = ns_patches.config.get_conf( subj_id = None,
                                           subj_types = "exp",
                                           excl_subj_ids = conf.ana.exclude_subj_ids
                                         )

    subj_ids = all_conf.all_subj.subj.keys()
    subj_ids.sort()

    n_subj = len( subj_ids )

    i_patches = range( conf.exp.n_mod_patches )

    # cull the patches with not enough nodes for all subjects
    i_valid_patches = np.setdiff1d( i_patches, conf.ana.exclude_patch_ids )

    data = np.empty( ( n_subj, 2 ) )
    data.fill( np.NAN )

    for ( i_subj, subj_id ) in enumerate( subj_ids ):

        subj_conf = ns_patches.config.get_conf( subj_id )
        subj_paths = ns_patches.paths.get_subj_paths( subj_conf )

        # this is patches x images x runs x type
        subj_data = np.load( subj_paths.ana.img_resp.full( ".npy" ) )

        # average over patches - note the indexing
        subj_data = np.mean( subj_data[ i_valid_patches, ... ], axis = 0 )
        # ... and images
        subj_data = np.mean( subj_data, axis = 0 )
        # ... and runs
        subj_data = np.mean( subj_data, axis = 0 )

        data[ i_subj, : ] = subj_data

    # save it in a human-readable form, since we dont need the precision
    np.savetxt( paths.coh_summ.full( ".txt" ), data, fmt = "%.5f" )

    # ... because we're going to bootstrap now

    coh_minus_incoh = data[ :, 0 ] - data[ :, 1 ]

    diff_mean = np.mean( coh_minus_incoh )

    diff_ci = scikits.bootstrap.ci( coh_minus_incoh )

    np.savetxt( paths.coh_diff_stats.full( ".txt" ),
                np.hstack( ( diff_mean, diff_ci ) ),
                fmt = "%.5f"
              )

