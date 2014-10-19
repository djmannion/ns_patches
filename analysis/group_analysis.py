
import numpy as np
import scipy.stats
import scikits.bootstrap

import ns_patches.config, ns_patches.paths


def loc_stats( conf, paths ):

    all_conf = ns_patches.config.get_conf( subj_id = None, subj_types = "loc" )

    subj_ids = all_conf.all_subj.subj.keys()
    subj_ids.sort()

    patch_k = np.empty( ( len( subj_ids ), all_conf.exp.n_mod_patches ) )
    patch_k.fill( np.NAN )

    for ( i_subj, subj_id ) in enumerate( subj_ids ):

        subj_conf = ns_patches.config.get_conf( subj_id )
        subj_paths = ns_patches.paths.get_subj_paths( subj_conf )

        subj_patch_k = np.loadtxt( subj_paths.loc.patch_id_count.full( ".txt" ) )

        patch_k[ i_subj, : ] = subj_patch_k

    assert np.sum( np.isnan( patch_k ) ) == 0

    np.savetxt( paths.patch_k.full( ".txt" ), patch_k, fmt = "%d" )


def coh_summ():

    conf = ns_patches.config.get_conf()
    group_paths = ns_patches.paths.get_group_paths()

    all_conf = ns_patches.config.get_conf(
        subj_id=None,
        subj_types="exp",
        excl_subj_ids=conf.ana.exclude_subj_ids
    )

    subj_ids = all_conf.all_subj.subj.keys()
    subj_ids.sort()

    n_subj = len(subj_ids)

    # cull the patches with not enough nodes for all subjects
    patch_ids = np.setdiff1d(
        conf.exp.mod_patches,
        conf.ana.exclude_patch_ids
    )

    patch_data = np.empty((n_subj, len(patch_ids), 2))
    patch_data.fill(np.NAN)

    for (i_subj, subj_id) in enumerate(subj_ids):

        subj_conf = ns_patches.config.get_conf( subj_id )
        subj_paths = ns_patches.paths.get_subj_paths( subj_conf )

        # nodes x (i, j, k, patch_id, patch_dist, all_t, coh, incoh)
        subj_data = np.loadtxt(subj_paths.coh_ana.comb.full(".txt"))

        i_act = subj_data[:, 5] > conf.ana.all_t_thresh

        for (i_patch, patch_id) in enumerate(patch_ids):

            # add 1 because they're stored as one-based
            i_patch_ok = subj_data[:, 3] == patch_id + 1

            assert np.sum(np.isnan(subj_data[i_patch_ok, -2:])) == 0

            patch_data[i_subj, i_patch, :] = np.mean(
                subj_data[i_patch_ok, -2:],
                axis=0
            )

    np.save(
        group_paths.coh_summ.full(".npy"),
        patch_data
    )

def depth_summ():

    conf = ns_patches.config.get_conf()
    group_paths = ns_patches.paths.get_group_paths()

    all_conf = ns_patches.config.get_conf(
        subj_id=None,
        subj_types="exp",
        excl_subj_ids=conf.ana.exclude_subj_ids
    )

    subj_ids = all_conf.all_subj.subj.keys()
    subj_ids.sort()

    n_subj = len(subj_ids)

    depth_data = np.empty((n_subj, len(conf.ana.bin_centres), 2))
    depth_data.fill(np.NAN)

    for (i_subj, subj_id) in enumerate(subj_ids):

        subj_conf = ns_patches.config.get_conf(subj_id)
        subj_paths = ns_patches.paths.get_subj_paths(subj_conf)

        for i_bin in xrange(len(conf.ana.bin_centres)):

            # this is i, j, k, patch, coh, non-coh
            subj_data = np.loadtxt(
                subj_paths.depth_ana.comb.full("_bin_{n:d}.txt".format(n=i_bin))
            )

            # average over all nodes and patches
            depth_data[i_subj, i_bin, :] = np.mean(subj_data[:, -2:], axis=0)

    assert np.sum(np.isnan(depth_data)) == 0

    np.save(group_paths.depth_summ.full(".npy"), depth_data)
