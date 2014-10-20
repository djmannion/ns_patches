
import numpy as np

import ns_patches.paths
import ns_patches.config


def patch_count_table(out_path):

    conf = ns_patches.config.get_conf()

    all_conf = ns_patches.config.get_conf(subj_id=None, subj_types="loc")

    subj_ids = all_conf.all_subj.subj.keys()
    subj_ids.sort()

    i_excluded_subj = subj_ids.index(
        conf.ana.exclude_subj_ids[0] + "_loc"
    )
    i_included_subj = np.setdiff1d(np.arange(len(subj_ids)), i_excluded_subj)

    group_paths = ns_patches.paths.get_group_paths()

    # this is subjects x modulated patches
    # it is sorted by subject ID
    # columns are indices into the modulated patches
    patch_k = np.loadtxt(group_paths.patch_k.full(".txt"))

    n_cols = patch_k.shape[1] + 4

    out = "\\begin{tabu}[c]{*{" + str(n_cols) + "}{c}}\n"

    header = (
        " & " +
        "\\multicolumn{" + str(patch_k.shape[0]) + "}{c}{Participant}" +
        " \\\\\n"
    )
    
    header += "\\cline{2-" + str(patch_k.shape[0] + 1) + "} \n"

    header += " & ".join(
        ["Patch"] +
        [str(x) for x in np.arange(1, patch_k.shape[0] + 1)] +
        ["Mean", "Min", "Min~(8/9)"]
    )

    out += header + " \\\\\n"

    out += "\\hline\n"

    for i_patch in xrange(patch_k.shape[1]):

        if (conf.exp.mod_patches[i_patch]) in conf.ana.exclude_patch_ids:
            text_colour = "red"
        else:
            text_colour = "black"

        # patch number
        entries = [str(conf.exp.mod_patches[i_patch] + 1)]

        for i_subj in xrange(patch_k.shape[0]):

            if i_subj == 5 or text_colour == "red":
                subj_colour = "red"
            else:
                subj_colour = "black"

            entries += [
                "\\color{" + subj_colour + "} " +
                str(int(patch_k[i_subj, i_patch]))
            ]

        entries += ["{n:.2f}".format(n=np.mean(patch_k[:, i_patch]))]
        entries += [str(int(np.min(patch_k[:, i_patch])))]
        entries += [str(int(np.min(patch_k[i_included_subj, i_patch])))]

        row = "\\rowfont{\\color{" + text_colour + "}} " + " & ".join(entries) + " \\\\ \n"

        out += row

#    out += "\\color{black} \n"

    out += "\\hline \n"

    subj_mean = np.mean(patch_k, axis=1)

    entries = (
        ["Mean"] +
        ["{n:.2f}".format(n=n) for n in subj_mean] +
        [""] * 3
    )

    out += " & ".join(entries) + " \\\\ \n"

    out += "\\end{tabu}"

    with open(out_path, "w") as out_file:
        out_file.writelines(out)
