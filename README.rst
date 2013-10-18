.. highlight:: tcsh

===================================================
Analysis code for natural scenes patches experiment
===================================================

Preprocessing
=============

The pre-processing pipelines are similar for the experiment and localiser sessions.
The subject ID has ``_loc`` appended to it for the localiser session.

Filesystem
~~~~~~~~~~

1. Make the subject's directory structure.
If the experiment session::

    mkdir -p ${SUBJ_ID}/{analysis,fmap,func/run{01,02,03,04,05,06,07,08},logs,reg}

If the localiser session::

    mkdir -p ${SUBJ_ID}/{analysis,fmap,func/run{01,02,03,04,05,06},reg}

2. If it is the experiment session, copy the subject's runtime logfiles to the ``logs`` directory.
If it is the localiser session, make a link to the standard timing files in the ``analysis`` directory::

    cd analysis
    ln -s ~/code/ns_patches/timing

3. Make symlinks named ``raw`` in each functional run directory that link to the location of its associated raw DICOM directory::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-such_and_such raw

4. Similarly, make symlinks named ``mag-raw`` and ``ph-raw`` in each fieldmap directory that link to the locations of the fieldmap acquisition::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-SEyada mag-raw
    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/PH-SEyada ph-raw

5. Make a local copy of the AFNI/SUMA base anatomical that we can use for alignment::

    3dcopy \
       {$SUBJECTS_DIR}/{$SUBJ_ID}/SUMA/{$SUBJ_ID}_SurfVol+orig \
       reg/{$SUBJ_ID}_ns_patches-anat+orig

The subject ID must be without ``_loc``, even for the localiser session, so it can match the reference anatomy.


Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files::

    ns_patches_preproc ${SUBJ_ID} convert

After execution, open up each NIFTI file and inspect for image quality and look at the summary image to see how much movement there was.


Slice-time correction
~~~~~~~~~~~~~~~~~~~~~

Adjusts for differences in slice acquisition timing::

    ns_patches_preproc ${SUBJ_ID} st_correct


Fieldmap preparation
~~~~~~~~~~~~~~~~~~~~

Prepares the fieldmap::

    ns_patches_preproc ${SUBJ_ID} fieldmap


Correction
~~~~~~~~~~
Applies a motion and distortion correction procedure::

    ns_patches_preproc ${SUBJ_ID} mc_unwarp

After execution, open up the summary NIFTI file to check that most of the motion has been removed.
To verify that the unwarping has worked correctly:

* Run ``fslview``.
* Load the original or corrected image from a given run.
* Add the magnitude image from the fieldmap as an overlay.
* Notice the geometric distortions in the functional data.
* Add the undistorted image as an overlay, and hide the uncorrected image.
* Toggle the visibility of the undistorted image, and verify that the geometry now aligns well with that of the fieldmap's magnitude image.


Anatomical registration
~~~~~~~~~~~~~~~~~~~~~~~

First, make a copy of the mean functional::

    cd reg
    3dcopy ../func/${SUBJ_ID}_ns_patches-mean.nii ${SUBJ_ID}_ns_patches-mean+orig

Again, if this is a localiser session, the subject ID can't have the ``_loc`` extension.

Now, we want to calculate some transformation parameters that will get the two images into rough register.
This will give the automated algorithm a good starting point.

* Start AFNI, from within the ``reg`` directory.
* Set the reference anatomical as the underlay.
* Position the crosshairs at a landmark on the brain. I like to use the most posterior portion of the occipital lobe, on the right side (in the image). Note down the three position values in the AFNI window (    in mm). Say they are ``[ 100, 50, 50 ]``.
* Then, change the underlay (or overlay, if you prefer) to the mean functional.
* Position the crosshairs at the same landmark as you used for the anatomical. The position might now be ``[ 20, 10, -20 ]``.
* Calculate ( reference anatomical positions - functional positions ), elementwise. In this example, that would give ``[ 80, 40, 70 ]``.
* Update the subject's configuration structure to include the estimate.

Then run::

    ns_patches_preproc ${SUBJ_ID} sess_reg


Surface projection
~~~~~~~~~~~~~~~~~~
Projects the functional images to a standardised cortical surface, averaging between the white and pial surfaces::

    ns_patches_preproc ${SUBJ_ID} vol_to_surf



Localiser analysis
==================

GLM
~~~

Estimate the GLM::

    ns_patches_loc_analysis sXXXX glm


Patch identification
~~~~~~~~~~~~~~~~~~~~

Before running, make a copy of the subject's visual localiser ROIs::

    cd analysis
    3dcopy \
      ${ROI_DIR}/${SUBJ_ID}_vis_loc_rois_${HEMI}-full.niml.dset \
      ${SUBJ_ID}_ns_patches-vis_loc_rois_${HEMI}-full.niml.dset

Work out the assignment of V1 nodes to patches::

    ns_patches_loc_analysis sXXXX patch_id


Patch cluster
~~~~~~~~~~~~~

Applies a small cluster threshold to get rid of salt-and-pepper noise::

    ns_patches_loc_analysis sXXXX patch_cluster


Patch count
~~~~~~~~~~~

Generates some stats (number of nodes) for each patch::

    ns_patches_loc_analysis sXXXX patch_count


Experiment analysis
===================

GLM
~~~

Before running, make a copy of the subject's visual localiser ROIs (we only run the GLM for V1 nodes)::

    cd analysis
    3dcopy \
      ${ROI_DIR}/${SUBJ_ID}_vis_loc_rois_${HEMI}-full.niml.dset \
      ${SUBJ_ID}_ns_patches-vis_loc_rois_${HEMI}-full.niml.dset

Estimate the GLM::

    ns_patches_exp_analysis sXXXX glm


Percent-signal-change conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run::

    ns_patches_exp_analysis sXXXX beta_to_psc


Patch dump
~~~~~~~~~~

This dumps the values for each image presentation for each node with an associated patch::

    ns_patches_exp_analysis sXXXX patch_dump


Image parcellate
~~~~~~~~~~~~~~~~

Distribute the patch information over image types::

    ns_patches_exp_analysis sXXXX image_resp

This results in ``${SUBJ_ID}_ns_patches-img_resp.npy``, which is a patch (26) x image (20) x run (8) x coherence (2) matrix.

