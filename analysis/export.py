"""Exports the timing so data can be analysed by others"""

import os

import numpy as np
import scipy.stats

import ns_patches.config
import ns_patches.paths
import ns_patches.exp.exp


def save_subj_timing():

    conf = ns_patches.config.get_conf(
        subj_types="exp"
    )
    paths = ns_patches.paths.get_exp_paths(conf)
    paths.img_db_info = os.path.join(
        "/home/damien/venv_study/ns_patches_analysis/",
        "code/ns_patches/ns_patches/img_info.txt"
    )

    log_dir = "/home/damien/venv_study/ns_patches_analysis/data/logs"
    timing_dir = "/home/damien/venv_study/ns_patches_analysis/data/exp_timing"

    img_info = ns_patches.exp.exp.load_img_info(conf, paths)

    valid_subj_ids = set(conf.all_subj.subj.keys()) - set(conf.ana.exclude_subj_ids)

    for subj_id in valid_subj_ids:

        run_data = [
            np.load(
                os.path.join(
                    log_dir,
                    "{s:s}_ns_patches-run_{n:02d}_log.npz".format(
                        s=subj_id,
                        n=run_num
                    )
                )
            )
            for run_num in xrange(1, conf.all_subj.subj[subj_id].n_runs + 1)
        ]

        # i_mod_patch is between 0 and 26, patch id is 26 values but between 0
        # and 32
        for (i_mod_patch, patch_id) in enumerate(conf.exp.mod_patches):

            for (i_img, img_id) in enumerate(img_info["id"]):

                for trial_type in ("coh", "noncoh"):

                    timing_path = os.path.join(
                        timing_dir,
                        "{s:s}_ns_patches-ap_{ap:d}_img_{im:d}_{t:s}.txt".format(
                            s=subj_id,
                            ap=patch_id,
                            im=img_id,
                            t=trial_type
                        )
                    )

                    timing_file = open(timing_path, "w")

                    for curr_run in run_data:

                        run_seq = curr_run["seq"]
                        img_trials = curr_run["img_trials"]

                        # pull out the info just for this patch
                        curr_info = img_trials[patch_id, :].astype("int")

                        # the ID of the 'coherent' image on each trial
                        coh_ids = scipy.stats.mode(
                            img_trials,
                            axis=0
                        )[0][0].astype("int")

                        for i_trial in xrange(len(curr_info)):

                            # current trial in this aperture shows current
                            # image
                            if curr_info[i_trial] == i_img:

                                if trial_type == "coh":
                                    if coh_ids[i_trial] == i_img:

                                        timing_file.write(
                                            "{n:.0f}\t".format(n=run_seq[i_trial])
                                        )

                                elif trial_type == "noncoh":
                                    if coh_ids[i_trial] != i_img:

                                        timing_file.write(
                                            "{n:.0f}\t".format(n=run_seq[i_trial])
                                        )
                        timing_file.write("\n")

                    timing_file.close()
