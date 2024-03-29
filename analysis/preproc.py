from __future__ import division

import os, os.path
import logging
import glob

import numpy as np

import fmri_tools.preproc, fmri_tools.utils
import runcmd


def convert( paths ):
    """Converts the functionals and fieldmaps from dicom to nifti"""

    logger = logging.getLogger( __name__ )
    logger.info( "Running conversion..." )

    # aggregate the dicom directories
    raw_paths = ( paths.func.raws +
                  [ paths.fmap.mag_raw, paths.fmap.ph_raw ]
                )

    # aggregate the image paths
    img_paths = ( paths.func.origs +
                  [ paths.fmap.mag, paths.fmap.ph ]
                )

    for ( raw_path, img_path ) in zip( raw_paths, img_paths ):

        fmri_tools.preproc.dcm_to_nii( dcm_path = raw_path.full(),
                                       nii_path = img_path.full( ".nii" )
                                     )


    # check that they are all unique
    are_unique = fmri_tools.utils.files_are_unique( [ img_path.full( ".nii" )
                                                      for img_path in img_paths
                                                    ]
                                                  )

    assert are_unique

    # files to go into the summary
    summ_paths = [ orig.full() for orig in paths.func.origs ]

    # make a summary image from the files
    fmri_tools.preproc.gen_sess_summ_img( epi_paths = summ_paths,
                                          out_path = paths.summ.orig.full(),
                                        )


def st_correct( conf, paths ):
    """Performs slice-time correction then shifts into the 'right' space"""

    logger = logging.getLogger( __name__ )
    logger.info( "Running slice-time correction..." )

    # first, need to work out the slice acquisition times. These will be the same
    # across all subjects, but might as well extract one for each subject.
    dcm_file = glob.glob( paths.func.raws[ 0 ].full( "MR*.dcm" ) )[ 0 ]

    acq_times = fmri_tools.utils.get_slice_acq_time_from_dcm( dcm_file )

    np.savetxt( paths.func.tpattern.full(), acq_times, fmt = "%.6f" )

    # now for the st correction
    for ( orig_path, st_path ) in zip( paths.func.origs, paths.func.sts ):

        fmri_tools.preproc.st_correct( orig_path = orig_path.full( ".nii" ),
                                       slice_path = st_path.full( ".nii" ),
                                       tpattern = "@" + paths.func.tpattern.full(),
                                       tr = conf.acq.tr_s,
                                       ignore_at_start = conf.loc.n_cull_vol
                                     )

        fmri_tools.preproc.reorient_img( in_path = st_path.full( ".nii" ),
                                         out_path = st_path.full( ".nii" ),
                                         reorient_str = conf.acq.ras
                                       )

    # files to go into the summary
    summ_paths = [ st.full() for st in paths.func.sts ]

    # make a summary image from the files
    fmri_tools.preproc.gen_sess_summ_img( epi_paths = summ_paths,
                                          out_path = paths.summ.st.full(),
                                        )


def fieldmaps( conf, paths ):
    """Prepare the fieldmaps"""

    logger = logging.getLogger( __name__ )
    logger.info( "Running fieldmap preparation..." )

    epi_run = conf.subj.mot_base - 1

    # set a corrected EPI to define the space to resample to
    ref_epi = paths.func.sts[ epi_run ].full( ".nii[0]" )

    # want to calculate a coarse brain mask from the epi
    mask_path = paths.fmap.mask.full( ".nii" )
    epi_path = paths.summ.st.full( ".nii[{n:d}]".format( n = epi_run ) )

    mask_cmd = [ "3dAutomask",
                 "-SI", "{n:d}".format( n = conf.subj.mask_SI ),
                 "-overwrite",
                 "-prefix", mask_path,
                 epi_path
               ]

    runcmd.run_cmd( " ".join( mask_cmd ) )

    fmri_tools.preproc.make_fieldmap( mag_path = paths.fmap.mag.full(),
                                      ph_path = paths.fmap.ph.full(),
                                      fmap_path = paths.fmap.fmap.full(),
                                      delta_te_ms = conf.acq.delta_te_ms,
                                      ref_img = ref_epi,
                                      recentre_ph = "mean",
                                      recentre_mask = mask_path,
                                      strip_params = [ "-surface_coil" ],
                                      strip_mag = False
                                    )


def mc_unwarp( conf, paths ):
    """Combined motion and distortion correction"""

    logger = logging.getLogger( __name__ )
    logger.info( "Running motion and distortion correction..." )

    st_paths = [ st_path.full() for st_path in paths.func.sts ]
    corr_paths = [ corr_path.full() for corr_path in paths.func.corrs ]
    uw_paths = [ uw_path.full() for uw_path in paths.func.uws ]

    i_run_base = conf.subj.mot_base - 1
    fugue_params = [ "--median" ]
    i_vol_base = conf.subj.vol_base

    # one subject's localisers were collected at 1.5mm resolution, and with different PE etc.
    if conf.subj.subj_id != "s1046_loc":
        uw_direction_fsl = conf.acq.ph_enc_dir
        uw_direction_afni = "LR"
        voxel_size = "1.0"
        dwell_ms = conf.acq.dwell_ms
    else:
        uw_direction_fsl = "z"
        uw_direction_afni = "SI"
        voxel_size = "1.5"
        dwell_ms = 0.71 / 2.0

    fmap_path = paths.fmap.fmap.full()

    fmri_tools.preproc.mc_unwarp( orig_paths = st_paths,
                                  corr_paths = corr_paths,
                                  uw_paths = uw_paths,
                                  fmap_path = fmap_path,
                                  dwell_ms = dwell_ms,
                                  uw_direction_fsl = uw_direction_fsl,
                                  uw_direction_afni = uw_direction_afni,
                                  voxel_size = voxel_size,
                                  i_vol_base = i_vol_base,
                                  i_run_base = i_run_base,
                                  fugue_params = fugue_params
                                )

    fmri_tools.preproc.gen_sess_summ_img( epi_paths = uw_paths,
                                          out_path = paths.summ.uw.full()
                                        )

    fmri_tools.preproc.mean_image( uw_paths, paths.summ.mean.full() )


def sess_reg( conf, paths ):
    """Coregisters the session anatomical with its mean functional"""

    logger = logging.getLogger( __name__ )

    logger.info( "Running registration..." )

    if conf.subj.extra_al_params is not None:
        extra_al_params = conf.subj.extra_al_params
    else:
        extra_al_params = None

    fmri_tools.preproc.img_reg( reg_dir = paths.reg.base.full(),
                                base_file = paths.reg.mean.file( "+orig" ),
                                mov_file = paths.reg.anat_ref.file( "+orig" ),
                                extra_al_params = extra_al_params,
                                epi_strip = "None" #"3dAutomask"
                              )



def vol_to_surf( conf, paths ):
    """Converts the functional volume-based images to SUMA surfaces."""

    logger = logging.getLogger( __name__ )
    logger.info( "Running volume to surface projection..." )

    start_dir = os.getcwd()

    for ( uw_file, surf_file, run_dir ) in zip( paths.func.uws,
                                                paths.func.surfs,
                                                paths.func.runs
                                              ):

        os.chdir( run_dir.full() )

        for hemi in [ "lh", "rh" ]:

            spec_file = paths.reg.spec.full( "_{hemi:s}.spec".format( hemi = hemi ) )

            # replace the subject ID with what FreeSurfer/SUMA considers the subject
            # ID to be
            spec_file = spec_file.replace( conf.subj.subj_id, conf.subj.fs_subj_id )

            surf_path = surf_file.full( "_{h:s}.niml.dset".format( h = hemi ) )

            surf_cmd = [ "3dVol2Surf",
                         "-spec", spec_file,
                         "-surf_A", "smoothwm",
                         "-surf_B", "pial",
                         "-map_func", "ave",
                         "-f_steps", "15",
                         "-f_index", "nodes",
                         "-sv", paths.reg.anat_reg.full( "+orig" ),
                         "-grid_parent", uw_file.full( ".nii" ),
                         "-out_niml", surf_path,
                         "-overwrite"
                       ]

            runcmd.run_cmd( " ".join( surf_cmd ) )

            # convert to full
            full_path = surf_file.full( "_{h:s}-full.niml.dset".format( h = hemi ) )

            node_str = "{n:d}".format( n = conf.subj.node_k[ hemi ] )
            fmri_tools.utils.sparse_to_full( in_dset = surf_path,
                                             out_dset = full_path,
                                             pad_node = node_str
                                           )

    os.chdir( start_dir )


def vol_to_surf_depth( conf, paths ):
    """Converts the functional volume-based images to SUMA surfaces."""

    logger = logging.getLogger( __name__ )
    logger.info( "Running volume to surface projection at different depths..." )

    start_dir = os.getcwd()

    for ( uw_file, surf_file, run_dir ) in zip( paths.func.uws,
                                                paths.func.surfs,
                                                paths.func.runs
                                              ):

        os.chdir( run_dir.full() )

        for hemi in [ "lh", "rh" ]:

            spec_file = paths.reg.spec.full( "_{hemi:s}.spec".format( hemi = hemi ) )

            # replace the subject ID with what FreeSurfer/SUMA considers the subject
            # ID to be
            spec_file = spec_file.replace( conf.subj.subj_id, conf.subj.fs_subj_id )

            for (i_bin, bin_centre) in enumerate(conf.ana.bin_centres):

                white_frac = bin_centre - (conf.ana.bin_width / 2.0)
                grey_frac = bin_centre + (conf.ana.bin_width / 2.0) - 1.0

                surf_path = surf_file.full(
                    "_bin_{i:d}_{h:s}.niml.dset".format( h = hemi, i = i_bin )
                )

                surf_cmd = [
                    "3dVol2Surf",
                    "-spec", spec_file,
                    "-surf_A", "smoothwm",
                    "-surf_B", "pial",
                    "-map_func", "ave",
                    "-f_steps", "15",
                    "-f_index", "nodes",
                    "-f_p1_fr", str(white_frac),
                    "-f_pn_fr", str(grey_frac),
                    "-sv", paths.reg.anat_reg.full( "+orig" ),
                    "-grid_parent", uw_file.full( ".nii" ),
                    "-out_niml", surf_path,
                    "-overwrite"
                ]

                runcmd.run_cmd( " ".join( surf_cmd ) )

                # convert to full
                full_path = surf_file.full(
                    "_bin_{i:d}_{h:s}-full.niml.dset".format(
                        h = hemi, i = i_bin
                    )
                )

                node_str = "{n:d}".format( n = conf.subj.node_k[ hemi ] )
                fmri_tools.utils.sparse_to_full( in_dset = surf_path,
                                                 out_dset = full_path,
                                                 pad_node = node_str
                                               )

    os.chdir( start_dir )

