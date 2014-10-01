
import os, os.path

import figutils
import matplotlib
figutils.set_defaults()

import matplotlib.pyplot as plt
plt.ioff()
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats

import fmri_tools.utils
import figutils

import ns_patches.config, ns_patches.paths, ns_patches.analysis.group_analysis
import ns_patches.analysis.analysis


def ecc_diff_fig(conf, paths):

    all_conf = ns_patches.config.get_conf(
            subj_id=None,
            subj_types="exp",
            excl_subj_ids=conf.ana.exclude_subj_ids
    )

    paths = ns_patches.paths.get_group_paths()

    coh_summ = ns_patches.analysis.group_analysis.patch_summ(all_conf, paths)

    patch_ids = np.setdiff1d(
        conf.exp.mod_patches,
        conf.ana.exclude_patch_ids
    )

    ring_ids = conf.stim.patches[patch_ids]["ring"]

    ring_diff = np.zeros(3)
    ring_err = np.zeros(3)

    for i_ring in xrange(3):

        ring_patches = ring_ids == i_ring

        ring_data = coh_summ[:, ring_patches]

        ring_diff[i_ring] = np.mean(ring_data)

        ring_mean = np.mean(ring_data, axis=1)

        ring_err[i_ring] = np.std(ring_mean, ddof=1) / np.sqrt(len(ring_mean))


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

    ring_ecc_degs = [1.8, 3.5, 6.1]

    for (i_ring, ring_ecc_deg) in enumerate(ring_ecc_degs):

        ax.plot(
            [ring_ecc_deg] * 2,
            [
                ring_diff[i_ring] - ring_err[i_ring],
                ring_diff[i_ring] + ring_err[i_ring]
            ],
            color="k",
            markersize=16
        )

        ax.plot([0,7], [0, 0], "k--")

        ax.hold(True)

        ax.plot(
            ring_ecc_deg,
            ring_diff[i_ring],
            markerfacecolor="k",
            marker="s",
            markeredgecolor="w",
            markersize=16,
            markeredgewidth=3
        )

    figutils.cleanup_fig(ax)

    ax.set_xlim(0, 7)
    ax.set_ylim(-0.02, 0.16)

#    ax.set_xlim(-1, 11.5)

#    ax.set_xticks([0.5]+range(3,coh_summ.shape[0] + 3))

#    ax.set_xticklabels(
#        ["Mean"] +
#        ["P{n:d}".format(n = n) for n in range(1, coh_summ.shape[0] + 1)]
#    )

    ax.set_xlabel("Patch eccenticity (deg visual angle)")
    ax.set_ylabel("Consistent - inconsistent (psc)")

    plt.show()


def coh_diff_fig(conf, paths):

    all_conf = ns_patches.config.get_conf(
            subj_id=None,
            subj_types="exp",
            excl_subj_ids=conf.ana.exclude_subj_ids
    )

    paths = ns_patches.paths.get_group_paths()

    coh_summ = ns_patches.analysis.group_analysis.patch_summ(all_conf, paths)

    coh_diff = np.mean(coh_summ, axis=1)

    std_err = np.std(coh_diff, ddof=1) / np.sqrt(len(coh_diff))

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



def plot_aperture_images(run_log_path, i_aperture=16, save_path=None):

    conf = ns_patches.config.get_conf()

    run_log = np.load(run_log_path)

    # this is (aperture x trials), where trials is between 80 and 88,
    # depending on how many null events there are in the pre-period
    img_trials = run_log["img_trials"]

    # need to find out which is the 'coherent' image on each trial
    (coh_ids, coh_counts) = scipy.stats.mode(img_trials, axis=0)

    assert np.all(
        coh_counts == (conf.stim.n_patches - conf.exp.n_incoh_patches)
    )

    coh_ids = coh_ids[0].astype("int")
    coh_ids = coh_ids[-conf.exp.n_trials:]

    # restrict it to the aperture
    ap_trials = img_trials[i_aperture, :]
    # and to the trial sequence
    ap_trials = ap_trials[-conf.exp.n_trials:]

    assert len(coh_ids) == len(ap_trials)

    run_mat = np.zeros((len(ap_trials), conf.exp.n_img))

    for (i_trial, (trial_img, trial_coh_img)) in enumerate(
        zip(ap_trials, coh_ids)
    ):

        if trial_img == trial_coh_img:
            trial_type = +1
        else:
            trial_type = -1

        run_mat[i_trial, trial_img] = trial_type

    fig = plt.figure(figsize=[3.3, 4], frameon=False)

    x_off = 0.13
    y_off = 0.15

    ax = plt.Axes(
            fig=fig,
            rect=[x_off, y_off, 0.97 - x_off, 0.97 - y_off],
            frameon=False
    )

    fig.add_axes(ax)

    ax.hold(True)

    ax.matshow(
        run_mat,
        aspect=0.3,
        interpolation="none",
        cmap="gray"
    )

    ax.set_xticks(np.arange(0, 21) - 0.5, minor=True)
    ax.set_yticks(np.arange(0, 81) - 0.5, minor=True)

    ax.set_xticklabels([])

    ax.tick_params(
            axis="both",
            direction="out",
            right="off",
            bottom="off",
            top="off"
    )

    ax.tick_params(
            axis="both",
            which="minor",
            left="off",
            right="off",
            top="off",
            bottom="off"
    )

    ax.grid(
        linestyle="-",
        linewidth=0.25,
        which="minor"
    )

    ax.set_xlabel("Source image")
    ax.set_ylabel("Run trial")

    patches = [
            mpatches.Patch(
                        facecolor=patch_colour,
                        edgecolor="black",
                        label=patch_label
                    )
            for (patch_colour, patch_label) in zip(
                        ("white", "black"),
                        ("Coherent", "Non-coherent")
                    )
    ]

    legend = plt.legend(
        handles=patches,
        ncol=2,
        loc=(0.2, -0.125)
    )

    legend.draw_frame(False)

    if save_path:
        plt.savefig(save_path)
