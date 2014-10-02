
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


def plot_cond_resp(save_path=None):

    conf = ns_patches.config.get_conf()
    group_paths = ns_patches.paths.get_group_paths()

    # this is (subj x aperture x condition [coh, non-coh])
    data = np.load(group_paths.coh_summ.full(".npy"))

    # average over apertures
    data = np.mean(data, axis=1)

    # do the stats now so can check normalisation
    (t, p) = scipy.stats.ttest_rel(data[:, 0], data[:, 1])

    # now normalise
    subj_mean = np.mean(data, axis=1)
    grand_mean = np.mean(data)

    data = (data - subj_mean[:, np.newaxis]) + grand_mean

    # this shouldn't affect the stats - check this is indeed true
    (norm_t, norm_p) = scipy.stats.ttest_rel(data[:, 0], data[:, 1])

    # floating point, so try and test safely
    np.testing.assert_almost_equal(t, norm_t)
    np.testing.assert_almost_equal(p, norm_p)

    # average over subjects
    cond_mean = np.mean(data, axis=0)
    # SEM, using sample stdev
    cond_sem = np.std(data, axis=0, ddof=1) / np.sqrt(data.shape[0])

    print "Condition means are: ", cond_mean
    print "Condition SEMs are: ", cond_sem
    print (
        "Difference stats are t(" +
        str(data.shape[0] - 1) +
        ") = " +
        str(t) +
        ", p = " +
        str(p)
    )

    # now can get on with the business of plotting!

    fig = plt.figure(figsize=[3.3, 2.5], frameon=False)

    x_off = 0.15
    y_off = 0.175
    x_max = 0.97
    y_max = 0.97
    y_lower_max = y_off + 0.15

    ax_base = plt.Axes(
        fig=fig,
        rect=[x_off, y_off, x_max - x_off, y_lower_max - y_off],
    )

    ax_plt = plt.Axes(
        fig=fig,
        rect=[x_off, y_lower_max, x_max - x_off, y_max - y_lower_max],
        sharex=ax_base,
    )

    fig.add_axes(ax_base)
    fig.add_axes(ax_plt)

    ax_plt.set_ylim([1.45, 1.65])
    ax_base.set_ylim([0,0.1])

    ax_base.set_xlim([-0.5, 1.5])

    ax_plt.plot(
        [0] * 2,
        [cond_mean[0] - cond_sem[0], cond_mean[0] + cond_sem[0]],
        "k"
    )

    ax_plt.plot(
        [1] * 2,
        [cond_mean[1] - cond_sem[1], cond_mean[1] + cond_sem[1]],
        "k"
    )

    ax_plt.scatter(
        [0, 1],
        cond_mean,
        facecolor=[0] * 3,
        edgecolor=[1] * 3,
        s=60,
        zorder=100
    )

    ax_plt.spines["bottom"].set_visible(False)
    ax_base.spines["top"].set_visible(False)

    ax_plt.tick_params(labeltop="off", labelbottom="off")
    ax_plt.tick_params(axis="x", bottom="off", top="off", right="off")
    ax_plt.tick_params(axis="y", right="off")
    ax_base.tick_params(axis="y", right="off")

    ax_base.spines["right"].set_visible(False)
    ax_plt.spines["right"].set_visible(False)
    ax_plt.spines["top"].set_visible(False)

    ax_base.xaxis.tick_bottom()

    ax_base.set_xlabel("Condition")
    ax_base.set_xticks([0, 1])
    ax_base.set_xticklabels(["Coherent", "Non-coherent"])

    ax_plt.set_ylabel("Response (normalised psc)", y=0.4)

    ax_base.spines["bottom"].set_position(("outward", 5))
    ax_base.spines["left"].set_position(("outward", 5))
    ax_plt.spines["left"].set_position(("outward", 5))

    ax_base.set_yticks([0])
    ax_base.tick_params(axis="y", length=0)

    kwargs = dict(transform=ax_base.transAxes, color='k', clip_on=False)

    ax_base.plot([-0.04, -0.01], [0.35, 0.45], "k", **kwargs)
    ax_base.plot([-0.04, -0.01], [0.45, 0.55], "k", **kwargs)

    if save_path:
        plt.savefig(save_path)

    plt.close(fig)


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

    plt.close(fig)
