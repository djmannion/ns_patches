import os, os.path

import fmri_tools.paths


def get_exp_paths( conf ):

    try:
        code_dir = os.environ[ "CODE_DIR" ]
    except KeyError:
        print "CODE_DIR not set in shell environment"
        raise

    exp = fmri_tools.paths.PathsHandler()

    exp.base_dir = os.path.join( code_dir, "ns_patches" )

    exp.timing_dir = os.path.join( exp.base_dir, "timing" )

    exp.img_db_info = os.path.join( exp.base_dir, "img_info.txt" )
    exp.img_db = os.path.join( exp.base_dir, "img_db" )

    exp.log_dir = os.path.join( exp.base_dir, "logs" )

    return exp


def get_subj_paths( conf ):
    """Get the path structure for a given subject"""

    base_dir = os.path.join( "/labs/olmanlab/Data7T/NatScenePatches/subj_data",
                             conf.subj.subj_id
                           )

    paths = fmri_tools.paths.get_default_paths( base_dir = base_dir,
                                                subj_id = conf.subj.subj_id,
                                                study_id = conf.exp.id,
                                                n_runs = conf.subj.n_runs
                                              )

    # add to func
    paths.func.tpattern = paths.func.base + ( conf.subj.subj_id + "_" +
                                              conf.exp.id + "-" +
                                              "tpattern.txt"
                                            )

    # add the st datatype
    paths.func.sts = [ orig + orig.file().replace( "orig", "st" )
                       for orig in paths.func.origs
                     ]

    paths.summ.st = ( paths.summ.orig +
                      paths.summ.orig.file().replace( "orig", "st" )
                    )

    paths.logs = _get_log_paths( conf, paths )

    if conf.subj.is_loc:
        paths.loc = _get_loc_paths( conf, paths )
        paths.reg = fmri_tools.paths.get_reg_paths( paths.reg.base.full(),
                                                    conf.subj.fs_subj_id,
                                                    conf.exp.id
                                                  )
    else:
        paths.ana = _get_ana_paths( conf, paths )

    paths.roi = _get_roi_paths( conf, paths )


    return paths


def _get_loc_paths( conf, paths ):
    """Get the paths for the localisers"""

    loc = fmri_tools.paths.PathsHandler()

    loc.base = paths.base / "analysis"

    subj_id = conf.subj.subj_id
    exp_id = conf.exp.id

    file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id,
                                                  exp_id = exp_id
                                                )

    loc.glm = loc.base + ( file_base + "glm" )
    loc.beta = loc.base + ( file_base + "beta" )
    loc.mask = loc.base + ( file_base + "mask" )

    loc.sig = loc.base + ( file_base + "sig" )
    loc.all_patch_id = loc.base + ( file_base + "all_patch_id" )
    loc.patch_id = loc.base + ( file_base + "patch_id" )
    loc.patch_id_count = loc.base + ( file_base + "patch_id_count" )
    loc.patch_id_thr = loc.base + ( file_base + "patch_id_thr" )
    loc.sig_sum = loc.base + ( file_base + "sig_sum" )
    loc.vl = loc.base + ( file_base + "vis_loc_rois" )

    loc.timing_base = ( loc.base /
                        "timing" +
                        "ns_patches-loc_timing_patch"
                      )

    return loc


def _get_ana_paths( conf, paths ):
    """Get the paths for the analysis"""

    ana = fmri_tools.paths.PathsHandler()

    ana.base = paths.base / "analysis"

    subj_id = conf.subj.subj_id
    exp_id = conf.exp.id

    file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id, exp_id = exp_id )

    ana.stim_times = ana.base + ( file_base + "stim_times" )

    ana.glm = ana.base + ( file_base + "glm" )
    ana.beta = ana.base + ( file_base + "beta" )

    ana.bltc = ana.base + ( file_base + "bltc" )
    ana.bl = ana.base + ( file_base + "bl" )
    ana.psc = ana.base + ( file_base + "psc" )

    ana.vl = ana.base + ( file_base + "vis_loc_rois" )
    ana.mask = ana.base + ( file_base + "mask" )

    ana.patch_resp = ana.base + ( file_base + "patch_resp" )

    ana.img_resp = ana.base + ( file_base + "img_resp" )

    ana.vec_resp = ana.base + ( file_base + "vec_resp" )

    ana.regress = ana.base + ( file_base + "regress_coef" )

    return ana


def _get_roi_paths( conf, paths ):
    """Get the paths for the ROI analysis"""

    roi = fmri_tools.paths.PathsHandler()

    roi.base = paths.base / "rois"

    subj_id = conf.subj.subj_id
    exp_id = conf.exp.id

    file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id, exp_id = exp_id )

    roi.rois = roi.base + ( file_base + "rois" )

    roi.psc = roi.base + ( file_base + "psc" )

    roi.psc_header = roi.base + ( file_base + "psc_header" )

    return roi


def _get_log_paths( conf, paths ):
    """Get the paths for the logfiles"""

    logs = fmri_tools.paths.PathsHandler()

    logs.base = paths.base / "logs"

    subj_id = conf.subj.subj_id
    exp_id = conf.exp.id

    file_base = "{subj_id:s}_{exp_id:s}-exp_log_".format( subj_id = subj_id, exp_id = exp_id )

    logs.logs = logs.base + ( file_base )

    logs.run_log_base = logs.base + ( subj_id + "_ns_patches-run_" )

    return logs


def get_group_paths( conf ):

    grp = fmri_tools.paths.PathsHandler()

    grp.base = fmri_tools.paths.Path( "/labs/olmanlab/Data7T/NatScenePatches/group_data" )

    file_base = "ns_patches"

    grp.patch_k = grp.base + ( file_base + "-patch_counts" )

    grp.coh_summ = grp.base + ( file_base + "-coh_summ" )
    grp.coh_diff_stats = grp.base + ( file_base + "-coh_diff_stats" )

    return grp




