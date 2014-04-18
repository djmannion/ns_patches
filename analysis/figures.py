
import os, os.path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

import fmri_tools.utils
import figutils

import ns_patches.config, ns_patches.paths, ns_patches.analysis.group_analysis


def coh_diff_fig(conf, paths):

    all_conf = ns_patches.config.get_conf(
            subj_id=None,
            subj_types="exp",
            excl_subj_ids=conf.ana.exclude_subj_ids
    )

    paths = ns_patches.paths.get_group_paths()

    coh_summ = ns_patches.analysis.group_analysis.patch_summ(all_conf, paths)

    coh_diff = np.mean(coh_summ, axis=1)

    std_err = np.std(coh_diff) / np.sqrt(len(coh_diff))

    ci = [np.mean(coh_diff), np.mean(coh_diff) - std_err, np.mean(coh_diff) + std_err]

    figutils.set_defaults()

    font_params = {
        "axes.labelsize": 24 * (1 / 1.25),
        "xtick.labelsize": 22 * (1 / 1.25),
        "ytick.labelsize": 22 * (1 / 1.25),
    }

    plt.rcParams.update(font_params)
    plt.ioff()

    fig = plt.figure()

    fig.set_size_inches(13, 9.3, forward = True)

    ax = plt.subplot(111)

    ms = 16

    ax.plot(
        range(3, coh_summ.shape[0] + 3),
        coh_diff,
        markerfacecolor="k",
        markeredgecolor="w",
        marker="o",
        linestyle="None",
        markersize=ms * 0.75
    )

    ax.hold(True)

    ax.plot(
        [0.5]*2,
        ci[1:],
        color="k",
        markersize=ms
    )

    ax.plot(
        [0.5],
        ci[0],
        markerfacecolor="k",
        marker="s",
        markeredgecolor="w",
        markersize=ms,
        markeredgewidth=3
    )

    ax.plot(
        [-1,11.5],
        [0, 0],
        color="k",
        linestyle="--"
    )

    figutils.cleanup_fig(ax)

    ax.set_xlim(-1, 11.5)

    ax.set_xticks([0.5]+range(3,coh_summ.shape[0] + 3))

    ax.set_xticklabels(
        ["Mean"] +
        ["P{n:d}".format(n = n) for n in range(1, coh_summ.shape[0] + 1)]
    )

    ax.set_xlabel("Participant")
    ax.set_ylabel("Consistent - inconsistent (psc)")

    plt.show()




def id_stats():

    all_conf = ns_patches.config.get_conf( subj_id = None, subj_types = "loc" )

    subj_ids = all_conf.all_subj.subj.keys()
    subj_ids.sort()

    figutils.set_defaults()

    fig = plt.figure()

    fig.set_size_inches( 7.08661, 10, forward = False )

    id_k = np.empty( ( len( subj_ids ), all_conf.exp.n_mod_patches ) )
    id_k.fill( np.NAN )

    for ( i_subj, subj_id ) in enumerate( subj_ids ):

        subj_conf = ns_patches.config.get_conf( subj_id )

        if not subj_conf.subj.is_loc:
            continue

        subj_paths = ns_patches.paths.get_subj_paths( subj_conf )

        k = np.loadtxt( subj_paths.loc.patch_id_count.full( ".txt" ) )

        id_k[ i_subj, : ] = k

        if np.any( k == 0 ):
            print "Subject " + subj_id + " has 0 node counts"

    print id_k

    assert np.sum( np.isnan( id_k ) ) == 0

    gs = gridspec.GridSpec( 1, 2 )

    for i_axis in xrange( 2 ):

        ax = plt.subplot( gs[ i_axis ] )

        mean_data = np.mean( id_k, axis = i_axis )
        err_data = np.std( id_k, axis = i_axis, ddof = 1 ) / np.sqrt( id_k.shape[ i_axis ] )

        ax.plot( range( 1, len( mean_data ) + 1 ), mean_data )

        figutils.cleanup_fig( ax )


#    plt.show()

    return id_k



def all_subj_scatter():

    all_conf = ns_patches.config.get_conf( None, True )

    subj_ids = all_conf.all_subj.subj.keys()

    _set_defaults()

    fig = plt.figure()

    fig.set_size_inches( 7.08661, 10, forward = True )

    gs = gridspec.GridSpec( 3, 3 )

    for ( i_subj, subj_id ) in enumerate( subj_ids ):

        conf = ns_patches.config.get_conf( subj_id )

        paths = ns_patches.paths.get_subj_paths( conf )

        ax = plt.subplot( gs[ i_subj ] )

        ax.hold( True )

        # patches x ( coh, non-coh )
        data = np.loadtxt( paths.ana.vec_resp.full( ".txt" ) )

        # slope, intercept
        coef = np.loadtxt( paths.ana.regress.full( ".txt" ) )

        ax.scatter( data[ :, 0 ],
                    data[ :, 1 ],
                    facecolors = "None",
                    edgecolors = "k"
                  )

        _cleanup_fig( ax )

        xlim = plt.xlim()
        ylim = plt.ylim()

        max = np.max( [ xlim, ylim ] )
        min = np.min( [ xlim, ylim ] )

        # axis line
        ax.plot( [ min, max ], [ 0, 0 ], color = [ 0.5 ] * 3 )
        ax.plot( [ 0, 0 ], [ min, max ], color = [ 0.5 ] * 3 )

        # unity line
        ax.plot( [ min, max ], [ min, max ], "b--" )

        fit_x = np.linspace( min, max, 100 )

        lin_fit = np.polyval( coef, fit_x )

        ax.plot( fit_x, lin_fit, "g" )

        ax.set_xlim( [ min, max ] )
        ax.set_ylim( [ min, max ] )

        ax.text( 0.4, 0.05,
                 ", ".join( [ subj_id,
                              "Slope: {s:.2f}".format( s = coef[ 0 ] ),
                              "Intercept: {s:.2f}".format( s = coef[ 1 ] )
                            ]
                          ),
                 transform = ax.transAxes,
                 fontsize = 8 / 1.25
               )

        if i_subj == len( subj_ids ) - 1:
            ax.set_xlabel( "Coherent response (psc)" )
            ax.set_ylabel( "Non-coherent response (psc)" )

    plt.subplots_adjust( left = 0.08,
                         bottom = 0.06,
                         right = 0.95,
                         top = 0.97,
                         wspace = 0.2,
                         hspace = 0.2
                       )

    plt.show()

