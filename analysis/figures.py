
import os
import functools

import figutils
import matplotlib
figutils.set_defaults()

import matplotlib.pyplot as plt
plt.ioff()
import matplotlib.patches as mpatches
import numpy as np
import scipy.stats
import svgutils.transform as sg

import fmri_tools.utils
import fmri_tools.stats

import ns_patches.config
import ns_patches.paths
import ns_patches.analysis.group_analysis
import ns_patches.analysis.analysis


def plot_depth(save_path=None):

    plot_depth_abs(save_path + "_abs.svg")
    plot_depth_diff(save_path + "_diff.svg")

    fig = sg.SVGFigure("8.38cm", "13.2cm")

    abs_fig = sg.fromfile(save_path + "_abs.svg")
    abs_plot = abs_fig.getroot()

    diff_fig = sg.fromfile(save_path + "_diff.svg")
    diff_plot = diff_fig.getroot()

    abs_plot.moveto(3.319, 11, scale=1.25)
    diff_plot.moveto(0, 247, scale=1.25)

    text = functools.partial(
        sg.TextElement,
        size=12,
        weight="bold",
        font="freesans"
    )

    A = text(11.6, 24, "A")
    B = text(11.6, 241, "B")

    fig.append([abs_plot, diff_plot])
    fig.append(A)
    fig.append(B)

    fig.save(save_path + ".svg")

    figutils.figutils.svg_to_pdf(
        svg_path=save_path + ".svg",
        pdf_path=save_path + ".pdf"
    )


def plot_depth_abs(save_path=None):

    conf = ns_patches.config.get_conf()
    group_paths = ns_patches.paths.get_group_paths()

    # this is subjects x depth x (coh, incoh)
    depth_data = np.load(group_paths.depth_summ.full(".npy"))

    subj_mean = np.mean(np.mean(depth_data, axis=-1), axis=-1)

    data_norm = (
        (depth_data - subj_mean[:, np.newaxis, np.newaxis]) +
        np.mean(depth_data)
    )

    data_mean = np.mean(data_norm, axis=0)
    data_sem = np.std(data_norm, axis=0, ddof=1) / np.sqrt(data_norm.shape[0])

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

    ax_plt.set_ylim([0.9, 2.1])
    ax_base.set_ylim([0, 0.1])

    symbols = ["s", "D"]
    labels = ["Coherent", "Non-coherent"]

    for (i_cond, flag) in enumerate([-1, +1]):

        ax_plt.plot(
            conf.ana.bin_centres + flag * 0.01,
            data_mean[:, i_cond],
            "k"
        )

        ax_plt.scatter(
            conf.ana.bin_centres + flag * 0.01,
            data_mean[:, i_cond],
            facecolor=[0] * 3,
            edgecolor=[1] * 3,
            s=60,
            zorder=100,
            marker=symbols[i_cond],
            label=labels[i_cond]
        )

        for (i_bin, bin_centre) in enumerate(conf.ana.bin_centres):

            ax_plt.plot(
                [bin_centre + flag * 0.01] * 2,
                [
                    data_mean[i_bin, i_cond] - data_sem[i_bin, i_cond],
                    data_mean[i_bin, i_cond] + data_sem[i_bin, i_cond]
                ],
                "k"
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

    ax_base.set_xlabel("Cortical depth (relative distance)")

    ax_base.set_xlim([-0.1, 1.1])

    ax_plt.set_ylabel("Response (normalised psc)", y=0.4)

    ax_base.spines["bottom"].set_position(("outward", 5))
    ax_base.spines["left"].set_position(("outward", 5))
    ax_plt.spines["left"].set_position(("outward", 5))

    ax_base.set_yticks([0])
    ax_base.tick_params(axis="y", length=0)

    ax_plt.set_xticks(conf.ana.bin_centres)

    xtick_labels = [
        "{n1:.1f} - {n2:.1f}".format(
            n1=bin_centre - conf.ana.bin_width / 2.0,
            n2=bin_centre + conf.ana.bin_width / 2.0
        )
        for bin_centre in conf.ana.bin_centres
    ]
    ax_plt.set_xticklabels(xtick_labels)

    kwargs = dict(transform=ax_base.transAxes, color='k', clip_on=False)

    ax_base.plot([-0.04, -0.01], [0.35, 0.45], "k", **kwargs)
    ax_base.plot([-0.04, -0.01], [0.45, 0.55], "k", **kwargs)

    leg = plt.legend(
        scatterpoints=1,
        loc="upper left"
    )

    leg.draw_frame(False)

    if save_path:
        plt.savefig(save_path)

    plt.close(fig)

    return leg


def plot_depth_diff(save_path=None):

    conf = ns_patches.config.get_conf()
    group_paths = ns_patches.paths.get_group_paths()

    depth_data = np.load(group_paths.depth_summ.full(".npy"))

    # coh - incoh
    depth_diff = depth_data[..., 0] - depth_data[..., 1]

    cwd = os.getcwd()

    stats = fmri_tools.stats.anova(depth_diff, "/home/damien/tmp", ["depth"])

    os.chdir(cwd)

    print stats

    # do stats on the difference
    for (i_bin, bin_centre) in enumerate(conf.ana.bin_centres):

        (t, p) = scipy.stats.ttest_1samp(depth_diff[:, i_bin], 0)

        print (
            "Bin centre at " + str(bin_centre) +
            ", " + str(t) + "; " + str(p)
        )

    trend_coeffs = [
        ("linear",  [-2, -1, +0, +1, +2]),
        ("quad",    [+2, -1, -2, -1, +2]),
        ("cubic",   [-1, +2, +0, -2, +1]),
        ("quartic", [+1, -4, +6, -4, +1])
    ]

    for (trend_name, trend_coeff) in trend_coeffs:

        (t, p) = scipy.stats.ttest_1samp(
            np.sum(depth_diff * trend_coeff, axis=1),
            0.0
        )

        print trend_name + " trend: " + str(t) + "; " + str(p)

    subj_mean = np.mean(depth_diff, axis=1)

    diff_norm = (depth_diff - subj_mean[:, np.newaxis]) + np.mean(depth_diff)
    diff_norm = depth_diff

    diff_mean = np.mean(diff_norm, axis=0)
    diff_sem = np.std(diff_norm, axis=0, ddof=1) / np.sqrt(diff_norm.shape[0])

    x_off = 0.16
    y_off = 0.175
    x_max = 0.97
    y_max = 0.97

    fig = plt.figure(figsize=[3.3, 2.5], frameon=False)

    ax_plt = plt.Axes(
        fig=fig,
        rect=[x_off, y_off, x_max - x_off, y_max - y_off]
    )

    fig.add_axes(ax_plt)

    ax_plt.plot(
        conf.ana.bin_centres,
        diff_mean
    )
    ax_plt.plot(
        [-0.05, 1.05],
        [0, 0],
        "--",
        color=[0.5] * 3
    )

    ax_plt.plot(
        conf.ana.bin_centres,
        diff_mean,
        "k"
    )

    ax_plt.scatter(
        conf.ana.bin_centres,
        diff_mean,
        facecolor=[0] * 3,
        edgecolor=[1] * 3,
        s=60,
        zorder=100
    )

    for (i_bin, bin_centre) in enumerate(conf.ana.bin_centres):

        ax_plt.plot(
            [bin_centre] * 2,
            [
                diff_mean[i_bin] - diff_sem[i_bin],
                diff_mean[i_bin] + diff_sem[i_bin]
            ],
            "k"
        )

    ax_plt.set_xlim([-0.1, 1.1])

    ax_plt.spines["bottom"].set_visible(True)

    ax_plt.tick_params(axis="x", top="off", right="off")
    ax_plt.tick_params(axis="y", right="off")

    ax_plt.spines["right"].set_visible(False)
    ax_plt.spines["top"].set_visible(False)

    ax_plt.set_xlabel("Cortical depth (relative distance)")
    ax_plt.set_ylabel("Coherent - non-coherent response (psc)")

    ax_plt.spines["bottom"].set_position(("outward", 5))
    ax_plt.spines["left"].set_position(("outward", 5))

    ax_plt.set_xticks(conf.ana.bin_centres)

    xtick_labels = [
        "{n1:.1f} - {n2:.1f}".format(
            n1=bin_centre - conf.ana.bin_width / 2.0,
            n2=bin_centre + conf.ana.bin_width / 2.0
        )
        for bin_centre in conf.ana.bin_centres
    ]
    ax_plt.set_xticklabels(xtick_labels)

    if save_path:
        plt.savefig(save_path)

    plt.close(fig)


def plot_dist(save_path=None):

    plot_dist_abs(save_path + "_abs.svg")
    plot_dist_diff(save_path + "_diff.svg")

    fig = sg.SVGFigure("8.38cm", "17.5cm")

    abs_fig = sg.fromfile(save_path + "_abs.svg")
    abs_plot = abs_fig.getroot()

    diff_fig = sg.fromfile(save_path + "_diff.svg")
    diff_plot = diff_fig.getroot()

    map_fig = sg.fromfile(
        os.path.join(
            "/home/damien/science/papers/current/ns_patches/figures",
            "dist/dist_map.svg"
        )
    )
    map_plot = map_fig.getroot()

    map_plot.moveto(41, 14, scale=1)
    abs_plot.moveto(0, 160.0, scale=1.25)
    diff_plot.moveto(0, 395.08, scale=1.25)

    text = functools.partial(
        sg.TextElement,
        size=12,
        weight="bold",
        font="freesans"
    )

    A = text(11.6, 13, "A")
    B = text(11.6, 162, "B")
    C = text(11.6, 390, "C")

    fig.append([abs_plot, diff_plot, map_plot])
    fig.append(A)
    fig.append(B)
    fig.append(C)

    fig.save(save_path + ".svg")

    figutils.figutils.svg_to_pdf(
        svg_path=save_path + ".svg",
        pdf_path=save_path + ".pdf"
    )


def plot_dist_diff(save_path=None):

    bin_start = 0.0
    bin_spacing = 2.0
    n_bins = 5

    bins = np.arange(0, n_bins + 1) * bin_spacing

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

    data = np.empty((n_subj, 2, n_bins + 1))
    data.fill(np.NAN)

    for (i_subj, subj_id) in enumerate(subj_ids):

        subj_data = np.loadtxt(
            os.path.join(
                group_paths.base.full(),
                subj_id + "_ns_patches-comb.txt"
            )
        )

        dist = subj_data[:, 4]

        for i_bin in xrange(n_bins):

            left_edge = i_bin * bin_spacing + bin_start
            right_edge = left_edge + bin_spacing

            in_bin = np.logical_and(
                dist >= left_edge,
                dist < right_edge
            )

            data[i_subj, :, i_bin] = np.mean(
                subj_data[in_bin, -2:],
                axis=0
            )

        in_bin = dist >= n_bins * bin_spacing

        data[i_subj, :, -1] = np.mean(
            subj_data[in_bin, -2:],
            axis=0
        )

    diff_data = data[:, 0, :] - data[:, 1, :]
    diff_mean = np.mean(diff_data, axis=0)
    diff_sem = np.std(diff_data, axis=0, ddof=1) / np.sqrt(diff_data.shape[0])

    cwd = os.getcwd()

    stats = fmri_tools.stats.anova(diff_data, "/home/damien/tmp", ["depth"])

    os.chdir(cwd)

    print stats

    trend_coeffs = [
        ("linear",  [-5, -3, -1, +1, +3, +5]),
        ("quad",    [+5, -1, -4, -4, -1, +5])
    ]

    for (trend_name, trend_coeff) in trend_coeffs:

        (t, p) = scipy.stats.ttest_1samp(
            np.sum(diff_data * trend_coeff, axis=1),
            0.0
        )

        print trend_name + " trend: " + str(t) + "; " + str(p)

    fig = plt.figure(figsize=[3.3, 2.5], frameon=False)

    x_off = 0.16
    y_off = 0.175
    x_max = 0.97
    y_max = 0.97

    ax_plt = plt.Axes(
        fig=fig,
        rect=[x_off, y_off, x_max - x_off, y_max - y_off]
    )

    fig.add_axes(ax_plt)

    ax_plt.set_ylim([-0.05, 0.2])
    ax_plt.set_xlim([-0.5, n_bins + 0.5])

    ax_plt.plot(
        [-0.5, 11.5],
        [0, 0],
        "--",
        color=[0.5] * 3
    )

    ax_plt.plot(
        np.arange(n_bins + 1),
        diff_mean,
        "k"
    )

    ax_plt.scatter(
        np.arange(n_bins + 1),
        diff_mean,
        facecolor=[0] * 3,
        edgecolor=[1] * 3,
        s=60,
        zorder=100
    )

    for i_bin in xrange(n_bins + 1):

        ax_plt.plot(
            [i_bin] * 2,
            [
                diff_mean[i_bin] - diff_sem[i_bin],
                diff_mean[i_bin] + diff_sem[i_bin]
            ],
            "k"
        )

    ax_plt.spines["bottom"].set_visible(True)

#    ax_plt.tick_params(labeltop="off", labelbottom="off")
    ax_plt.tick_params(axis="x", top="off", right="off")
    ax_plt.tick_params(axis="y", right="off")

    ax_plt.spines["right"].set_visible(False)
    ax_plt.spines["top"].set_visible(False)

    ax_plt.set_xlabel("Distance from aperture centre (mm)")
    ax_plt.set_xticks(range(n_bins + 1))

    xtick_labels = [
        "{n1:.0f} - {n2:.0f}".format(n1=bins[n1], n2=bins[n1] + bin_spacing)
        for n1 in np.arange(n_bins + 1)
    ]
    xtick_labels[-1] = "10+"

    ax_plt.set_xticklabels(xtick_labels)

    ax_plt.set_ylabel("Coherent - non-coherent response (psc)")

    ax_plt.spines["bottom"].set_position(("outward", 5))
    ax_plt.spines["left"].set_position(("outward", 5))

    if save_path:
        plt.savefig(save_path)

    plt.close(fig)


def plot_dist_abs(save_path=None):

    bin_start = 0.0
    bin_spacing = 2.0
    n_bins = 5

    bins = np.arange(0, n_bins + 1) * bin_spacing

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

    data = np.empty((n_subj, 2, n_bins + 1))
    data.fill(np.NAN)

    for (i_subj, subj_id) in enumerate(subj_ids):

        subj_data = np.loadtxt(
            os.path.join(
                group_paths.base.full(),
                subj_id + "_ns_patches-comb.txt"
            )
        )

        dist = subj_data[:, 4]

        for i_bin in xrange(n_bins):

            left_edge = i_bin * bin_spacing + bin_start
            right_edge = left_edge + bin_spacing

            in_bin = np.logical_and(
                dist >= left_edge,
                dist < right_edge
            )

            data[i_subj, :, i_bin] = np.mean(
                subj_data[in_bin, -2:],
                axis=0
            )

        in_bin = dist >= n_bins * bin_spacing

        data[i_subj, :, -1] = np.mean(
            subj_data[in_bin, -2:],
            axis=0
        )

    subj_mean = np.mean(np.mean(data, axis=-1), axis=-1)
    grand_mean = np.mean(data)

    data = (data - subj_mean[:, np.newaxis, np.newaxis]) + grand_mean

    data_mean = np.mean(data, axis=0)
    data_sem = np.std(data, axis=0, ddof=1) / np.sqrt(data.shape[0])

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

    ax_plt.set_ylim([0.7, 2.1])
    ax_base.set_ylim([0, 0.1])

    ax_base.set_xlim([-0.5, n_bins + 0.5])

    symbols = ["s", "D"]
    labels = ["Coherent", "Non-coherent"]

    for (i_cond, flag) in enumerate([-1, +1]):

        ax_plt.plot(
            np.arange(n_bins + 1) + flag * 0.1,
            np.mean(data[:, i_cond, :], axis=0),
            "k"
        )

        ax_plt.scatter(
            np.arange(n_bins + 1) + flag * 0.1,
            np.mean(data[:, i_cond, :], axis=0),
            facecolor=[0] * 3,
            edgecolor=[1] * 3,
            s=60,
            zorder=100,
            marker=symbols[i_cond],
            label=labels[i_cond]
        )

        for i_bin in xrange(n_bins + 1):

            ax_plt.plot(
                [i_bin + flag * 0.1] * 2,
                [
                    data_mean[i_cond, i_bin] - data_sem[i_cond, i_bin],
                    data_mean[i_cond, i_bin] + data_sem[i_cond, i_bin]
                ],
                "k"
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

    ax_base.set_xlabel("Distance from aperture centre (mm)")
    ax_base.set_xticks(range(n_bins + 1))
    ax_base.set_xticklabels(["Inner", "Middle", "Outer"])

    xtick_labels = [
        "{n1:.0f} - {n2:.0f}".format(n1=bins[n1], n2=bins[n1] + bin_spacing)
        for n1 in np.arange(n_bins + 1)
    ]
    xtick_labels[-1] = "10+"

    ax_base.set_xticklabels(xtick_labels)

    ax_plt.set_ylabel("Response (normalised psc)", y=0.4)

    ax_base.spines["bottom"].set_position(("outward", 5))
    ax_base.spines["left"].set_position(("outward", 5))
    ax_plt.spines["left"].set_position(("outward", 5))

    ax_base.set_yticks([0])
    ax_base.tick_params(axis="y", length=0)

    kwargs = dict(transform=ax_base.transAxes, color='k', clip_on=False)

    ax_base.plot([-0.04, -0.01], [0.35, 0.45], "k", **kwargs)
    ax_base.plot([-0.04, -0.01], [0.45, 0.55], "k", **kwargs)

    leg = plt.legend(
        scatterpoints=1,
        loc="upper right"
    )

    leg.draw_frame(False)

    if save_path:
        plt.savefig(save_path)

    plt.close(fig)


def plot_cond_resp_by_ecc(save_path=None):

    conf = ns_patches.config.get_conf()
    group_paths = ns_patches.paths.get_group_paths()

    # this is (subj x aperture x condition [coh, non-coh])
    data = np.load(group_paths.coh_summ.full(".npy"))

    # cull the patches with not enough nodes for all subjects
    patch_ids = np.setdiff1d(
        conf.exp.mod_patches,
        conf.ana.exclude_patch_ids
    )

    assert len(patch_ids) == data.shape[1]

    rings = [[], [], []]

    for (i_patch, patch_id) in enumerate(patch_ids):

        i_stim_patch = np.where(
            conf.stim.patches["id"] == patch_id
        )[0][0]

        i_ring = conf.stim.patches[i_stim_patch]["ring"]

        rings[i_ring].append(i_patch)

    # subj x coh x ring
    ecc_data = np.empty((data.shape[0], 2, 3))
    ecc_data.fill(np.NAN)

    for i_ring in xrange(3):
        ecc_data[..., i_ring] = np.mean(data[:, rings[i_ring], :], axis=1)

    subj_mean = np.mean(np.mean(ecc_data, axis=-1), axis=-1)
    grand_mean = np.mean(ecc_data)

    norm_data = (ecc_data - subj_mean[:, np.newaxis, np.newaxis]) + grand_mean

    mean_data = np.mean(norm_data, axis=0)
    sem_data = np.std(norm_data, axis=0, ddof=1) / np.sqrt(data.shape[0])

    cwd = os.getcwd()

    s = fmri_tools.stats.anova(
        ecc_data,
        "/home/damien/tmp",
        ["coh", "ring"]
    )

    print s

    s2 = fmri_tools.stats.anova(
        ecc_data[:, 0, :] - ecc_data[:, 1, :],
        "/home/damien/tmp",
        ["ring"]
    )

    print s2

    os.chdir(cwd)

    # do the stats
    for i_ring in xrange(3):

        print "Ring: " + str(i_ring + 1)
        print "\tMean: ", mean_data[:, i_ring]
        print "\tSEM: ", sem_data[:, i_ring]

        (t, p) = scipy.stats.ttest_rel(
            norm_data[:, 0, i_ring],
            norm_data[:, 1, i_ring]
        )

        print "\tt: ", t
        print "\tp: ", p

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

    ax_plt.set_ylim([1.0, 2.05])
    ax_base.set_ylim([0, 0.1])

    ax_base.set_xlim([-0.5, 2.5])

    symbols = ["s", "D"]
    labels = ["Coherent", "Non-coherent"]

    for (i_cond, flag) in enumerate([-1, +1]):

        ax_plt.scatter(
            np.arange(3) + flag * 0.1,
            mean_data[i_cond, :],
            facecolor=[0] * 3,
            edgecolor=[1] * 3,
            s=60,
            zorder=100,
            marker=symbols[i_cond],
            label=labels[i_cond]
        )

        for i_x in xrange(3):

            ax_plt.plot(
                [i_x + flag * 0.1] * 2,
                [
                    mean_data[i_cond, i_x] - sem_data[i_cond, i_x],
                    mean_data[i_cond, i_x] + sem_data[i_cond, i_x]
                ],
                "k"
            )

    for i_x in [1, 2]:
        ax_plt.plot([i_x - 0.1, i_x + 0.1], [1.975] * 2, "k")
        ax_plt.text(i_x, 1.98, "*", horizontalalignment="center")

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

    ax_base.set_xlabel("Aperture eccentricity")
    ax_base.set_xticks(range(3))
    ax_base.set_xticklabels(["Inner", "Middle", "Outer"])

    ax_plt.set_ylabel("Response (normalised psc)", y=0.4)

    ax_base.spines["bottom"].set_position(("outward", 5))
    ax_base.spines["left"].set_position(("outward", 5))
    ax_plt.spines["left"].set_position(("outward", 5))

    ax_base.set_yticks([0])
    ax_base.tick_params(axis="y", length=0)

    kwargs = dict(transform=ax_base.transAxes, color='k', clip_on=False)

    ax_base.plot([-0.04, -0.01], [0.35, 0.45], "k", **kwargs)
    ax_base.plot([-0.04, -0.01], [0.45, 0.55], "k", **kwargs)

    leg = plt.legend(
        scatterpoints=1,
        loc="upper left"
    )

    leg.draw_frame(False)

    if save_path:
        plt.savefig(save_path + ".svg")

        figutils.figutils.svg_to_pdf(
            svg_path=save_path + ".svg",
            pdf_path=save_path + ".pdf"
        )

    plt.close(fig)


def plot_cond_resp(save_path=None):

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
    ax_base.set_ylim([0, 0.1])

    ax_base.set_xlim([-0.5, 1.5])

    for i_cond in xrange(2):
        ax_plt.plot(
            [i_cond] * 2,
            [
                cond_mean[i_cond] - cond_sem[i_cond],
                cond_mean[i_cond] + cond_sem[i_cond]
            ],
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
        plt.savefig(save_path + ".svg")

        figutils.figutils.svg_to_pdf(
            svg_path=save_path + ".svg",
            pdf_path=save_path + ".pdf"
        )

    plt.close(fig)


def plot_aperture_images(run_log_path=None, i_aperture=16, save_path=None):

    if run_log_path is None:
        run_log_path = os.path.join(
            "/home/damien/venv_study/ns_patches_analysis/data/",
            "s1021_ns_patches-run_01_log.npz"
        )

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

    fig = plt.figure(
        figsize=[3.3, 4],
        frameon=False
    )

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
            label=patch_label,
            rasterized=False,
            aa=False
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

        plt.savefig(save_path + ".svg")

        figutils.figutils.svg_to_pdf(
            svg_path=save_path + ".svg",
            pdf_path=save_path + ".pdf"
        )

    plt.close(fig)
