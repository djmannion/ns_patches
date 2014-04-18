"""Performs the analysis for a single subject"""

import os
import os.path
import logging

import numpy as np
import scipy.stats

import fmri_tools.utils
import runcmd

import ns_patches.config
import ns_patches.paths


def _load_run_log(run_num, paths):
    "Loads the log file for a particular subject and run"

    run_ext = "{r:02d}_log.npz".format(r=run_num)
    run_log = np.load(paths.logs.run_log_base.full(run_ext))

    return run_log


def _get_coh_timing(conf, paths, patch_id):
    """Determine the timing for the coherent and incoherent events for a given
    patch"""

    coh_timing = []
    incoh_timing = []

    for run_num in xrange(1, conf.subj.n_runs + 1):

        coh_run_timing = []
        incoh_run_timing = []

        run_log = _load_run_log(run_num, paths)

        # n_trials long
        run_seq = run_log["seq"]

        # patch IDs x trials - each shows the image ID
        img_trials = run_log["img_trials"]

        # the patch IDs in `img_trials` are for ALL patches - just check to
        # make sure
        assert img_trials.shape[0] == conf.stim.n_patches

        # the ID of the 'coherent' image on each trial
        coh_ids = scipy.stats.mode(
            img_trials,
            axis=0
        )[0][0].astype("int")

        # just check we've pulled out what we think we have
        assert len(coh_ids) == len(run_seq)

        # extract the image IDs for the desired patch
        patch_img_ids = img_trials[patch_id, :].astype("int")

        # loop through each trial, pickout out the coherent image id, the patch
        # id, and the trial onset
        for (
            coh_img_id,
            patch_img_id,
            trial_seq
        ) in zip(
            coh_ids,
            patch_img_ids,
            run_seq
        ):

            # which timing list it gets assigned to depends on whether it is
            # the coherent image or not
            if coh_img_id == patch_img_id:
                coh_run_timing.append(trial_seq)
            else:
                incoh_run_timing.append(trial_seq)

        coh_timing.append(coh_run_timing)
        incoh_timing.append(incoh_run_timing)

    return (coh_timing, incoh_timing)


def coh_glm(conf, paths):
    "Run the GLM(s) for a given subject"

    # we only want to run GLMs for the patches with acceptable node counts, so
    # remove those identified as lacking
    patch_ids = np.setdiff1d(
        conf.exp.mod_patches,
        conf.ana.exclude_patch_ids
    )

    # run the GLM for each patch
    for patch_id in patch_ids:
        _run_coh_glm(conf, paths, patch_id)

    vf_lookup = {"lh": "R", "rh": "L"}

    # combine all into one
    for hemi in ["lh", "rh"]:

        comb_cmd = [
            "3dMean",
            "-non_zero",
            "-sum",
            "-prefix", paths.coh_ana.comb.full("_" + hemi + "-full.niml.dset"),
            "-overwrite"
        ]

        for patch_id in patch_ids:

            if conf.stim.patches[patch_id]["vf"] == vf_lookup[hemi]:

                comb_cmd.append(
                    paths.coh_ana.psc.file(
                        "-patch_{n:d}".format(n=patch_id) +
                        "_" + hemi +
                        "-full.niml.dset"
                    )
                )

        runcmd.run_cmd(" ".join(comb_cmd))


def _run_coh_glm(conf, paths, patch_id):
    "Run the coh/incoh GLM for a given patch"

    os.chdir(paths.coh_ana.base.full())

    # contralateral organisation
    if conf.stim.patches[patch_id]["vf"] == "L":
        hemi = "rh"
    else:
        hemi = "lh"

    hemi_ext = "_" + hemi

    # [coh, incoh]
    timings = _get_coh_timing(conf, paths, patch_id)

    # write patch timings
    for (cond, cond_name) in zip(timings, ["coh", "incoh"]):

        cond_path = paths.coh_ana.stim_times.full(
            "-patch_{n:d}_{c:s}.txt".format(n=patch_id, c=cond_name)
        )

        with open(cond_path, "w") as cond_file:

            for run_cond in cond:
                cond_file.write("\t".join(
                    [
                        "{t:d}".format(t=x) for x in run_cond
                    ]
                ))
                cond_file.write("\n")

    # write out the mask for this patch
    loc_id = conf.subj.subj_id + "_loc"
    loc_conf = ns_patches.config.get_conf(loc_id)
    loc_paths = ns_patches.paths.get_subj_paths(loc_conf)

    id_path = loc_paths.loc.patch_id_thr.full(
        "_{h:s}-full_Clustered_e1_a{n:.01f}.niml.dset".format(
            h=hemi,
            n=conf.loc.area_thr
        )
    )

    mask_path = paths.coh_ana.mask.full(
        "-patch_{n:d}".format(n=patch_id) +
        hemi_ext +
        "-full.niml.dset"
    )

    # the patch ID has a +1 in the below because they are stored in the niml as
    # 1-based
    mask_cmd = [
        "3dcalc",
        "-a", id_path,
        "-expr", "equals(a,{x:d})".format(x=patch_id + 1),
        "-prefix", mask_path,
        "-overwrite"
    ]

    runcmd.run_cmd(" ".join(mask_cmd))

    # right-o, ready for the GLM
    censor_vols = conf.exp.n_censor_vols - 1
    censor_str = "*:0-{v:.0f}".format(v=censor_vols)

    model_str = "SPMG1({d:.0f})".format(d=conf.exp.img_on_s)

    glm_cmd = [
        "3dDeconvolve",
        "-input"
    ]

    surf_paths = [
        surf_path.full(hemi_ext + "-full.niml.dset")
        for surf_path in paths.func.surfs
    ]

    glm_cmd.extend(surf_paths)

    glm_cmd.extend(
        [
            "-force_TR", "{tr:.3f}".format(tr=conf.acq.tr_s),
            "-polort", "a",  # auto baseline degree
            "-local_times",
            "-mask", mask_path,
            "-CENSORTR", censor_str,
            "-xjpeg", "exp_design_patch_{x:d}.png".format(x=patch_id),
            "-x1D", "exp_design_patch_{x:d}".format(x=patch_id),
            "-overwrite",
            "-x1D_stop",  # want to use REML, so don't bother running
            "-num_stimts", "2"
        ]
    )

    for (i_cond, cond_name) in enumerate(["coh", "incoh"]):

        glm_cmd.extend(
            [
                "-stim_label",
                "{x:d}".format(x=i_cond + 1), cond_name
            ]
        )

        glm_cmd.extend(
            [
                "-stim_times",
                "{x:d}".format(x=i_cond + 1),
                paths.coh_ana.stim_times.full(
                    "-patch_{n:d}_{c:s}.txt".format(n=patch_id, c=cond_name)
                ),
                model_str
            ]
        )

    # all > blank contrast
    con_str = "SYM: +coh +incoh"
    glm_cmd.extend(
        [
            "-gltsym", "'" + con_str + "'",
            "-glt_label", "1", "all_gt_bl"
        ]
    )

    runcmd.run_cmd(" ".join(glm_cmd))

    os.remove("Decon.REML_cmd")

    beta_file = paths.coh_ana.beta.full(
        "-patch_{n:d}".format(n=patch_id) +
        hemi_ext +
        "-full.niml.dset"
    )
    buck_file = paths.coh_ana.glm.full(
        "-patch_{n:d}".format(n=patch_id) +
        hemi_ext +
        "-full.niml.dset"
    )

    reml_cmd = [
        "3dREMLfit",
        "-matrix", "exp_design_patch_{x:d}.xmat.1D".format(x=patch_id),
        "-mask", mask_path,
        "-Rbeta", beta_file,
        "-tout",
        "-Rbuck", buck_file,
        "-overwrite",
        "-input"
    ]

    reml_cmd.append("'" + " ".join(surf_paths) + "'")

    # run the proper GLM
    runcmd.run_cmd(" ".join(reml_cmd))

    # now to convert to PSC, while we're here
    design_path = "exp_design_patch_{x:d}.xmat.1D".format(x=patch_id)

    # to write
    ext = "-patch_{n:d}".format(n=patch_id) + hemi_ext + "-full.niml.dset"
    bltc_path = paths.coh_ana.bltc.file(ext)
    bl_path = paths.coh_ana.bl.file(ext)
    psc_path = paths.coh_ana.psc.file(ext)

    # 4 orthogonal polynomial regressors per run
    n_nuisance = conf.subj.n_runs * 4

    # checked via '-verb'
    beta_bricks = "[{n:d}..$]".format(n=n_nuisance)

    fmri_tools.utils.beta_to_psc(
        beta_file,
        beta_bricks,
        design_path,
        bltc_path,
        bl_path,
        psc_path,
    )


def _patch_cent_dist(conf, paths, patch_id):
    "Write the distance of each node in a patch to its centre"

    os.chdir(paths.coh_ana.base.full())

    # contralateral organisation
    if conf.stim.patches[patch_id]["vf"] == "L":
        hemi = "rh"
    else:
        hemi = "lh"

    hemi_ext = "_" + hemi

    mask_path = paths.coh_ana.mask.full(
        "-patch_{n:d}".format(n=patch_id) +
        hemi_ext +
        "-full.niml.dset"
    )

    centre_node = fmri_tools.utils.get_centre_node(
        surf_dset=mask_path,
        spec_path=paths.reg.spec.full(hemi_ext + ".spec"),
        calc_surf="pial"
    )

    fmri_tools.utils.write_dist_to_centre(
        centre_node=centre_node,
        in_dset=mask_path,
        spec_path=paths.reg.spec.full(hemi_ext + ".spec"),
        dist_dset=paths.coh_ana.patch_dist.file(
            "-patch_{n:d}".format(n=patch_id) +
            hemi_ext
        ),
        pad_to=str(conf.subj.node_k[hemi]),
        inc_centre_node=True,
        calc_surf="pial"
    )


def centre_distances(conf, paths):
    "Write the distance of each node to its patch centre"

    for patch_id in conf.ana.valid_patch_ids:
        _patch_cent_dist(conf, paths, patch_id)

    os.chdir(paths.coh_ana.base.full())

    # now to combine
    vf_lookup = {"lh": "R", "rh": "L"}

    # combine all into one
    for hemi in ["lh", "rh"]:

        comb_cmd = [
            "3dMean",
            "-non_zero",
            "-sum",
            "-prefix", paths.coh_ana.patch_dist.full("_" + hemi + "-full.niml.dset"),
            "-overwrite"
        ]

        for patch_id in conf.ana.valid_patch_ids:

            if conf.stim.patches[patch_id]["vf"] == vf_lookup[hemi]:

                comb_cmd.append(
                    paths.coh_ana.patch_dist.file(
                        "-patch_{n:d}".format(n=patch_id) +
                        "_" + hemi +
                        ".niml.dset"
                    )
                )

        runcmd.run_cmd(" ".join(comb_cmd))



def data_dump(conf, paths):
    "Dump all the relevant data to a single text file"

    # ingredients:
    #   -node coordinates
    #   -patch IDp
    #   -patch centre dist
    #   -all > blank from GLM
    #   -coherent PSC
    #   -incoherent PSC

    patch_ids = np.setdiff1d(
        conf.exp.mod_patches,
        conf.ana.exclude_patch_ids
    )

    vf_lookup = {"lh": "R", "rh": "L"}

    # open it this way so that we can write to it twice and it will append
    dump_handle = open(paths.coh_ana.comb.full(".txt"), "w")

    os.chdir(paths.coh_ana.base.full())

    for hemi in ["lh", "rh"]:

        hemi_ext = "_" + hemi

        glm_comb_path = paths.coh_ana.glm_comb.full(
            hemi_ext + "-full.niml.dset"
        )

        # first, want to combine all the GLM data together
        comb_cmd = [
            "3dMean",
            "-non_zero",
            "-sum",
            "-prefix", glm_comb_path,
            "-overwrite"
        ]

        for patch_id in patch_ids:

            if conf.stim.patches[patch_id]["vf"] == vf_lookup[hemi]:

                buck_file = paths.coh_ana.glm.full(
                    "-patch_{n:d}".format(n=patch_id) +
                    hemi_ext +
                    "-full.niml.dset" +
                    "[6]"
                )

                comb_cmd.append(buck_file)

        runcmd.run_cmd(" ".join(comb_cmd))

        psc_comb_path = paths.coh_ana.comb.full(hemi_ext + "-full.niml.dset")

        # now we need to know the patch IDs
        loc_id = conf.subj.subj_id + "_loc"
        loc_conf = ns_patches.config.get_conf(loc_id)
        loc_paths = ns_patches.paths.get_subj_paths(loc_conf)

        id_path = loc_paths.loc.patch_id_thr.full(
            "_{h:s}-full_Clustered_e1_a{n:.01f}.niml.dset".format(
                h=hemi,
                n=conf.loc.area_thr
            )
        )

        dump_path = paths.coh_ana.comb.full(hemi_ext + ".txt")

        if os.path.exists(dump_path):
            os.remove(dump_path)

        # now we have all our ingredients, we can write out the text file
        dump_cmd = [
            "3dmaskdump",
            "-o", dump_path,
            "-mask", psc_comb_path,
            id_path,
            paths.coh_ana.patch_dist.full("_" + hemi + "-full.niml.dset"),
            glm_comb_path,
            psc_comb_path
        ]

        runcmd.run_cmd(" ".join(dump_cmd))

        # concatenate across hemis
        np.savetxt(dump_handle, np.loadtxt(dump_path))

    dump_handle.close()
