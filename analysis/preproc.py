"""
Set of routines to pre-process the fMRI data for the diamond
 fMRI experiment.
"""

from __future__ import division

import os, os.path
import logging
import glob

import numpy as np
import scipy.io

import nipy.io.api
import nipy.algorithms.resample

import fmri_tools.preproc, fmri_tools.utils


def convert( conf, paths ):
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
	assert( fmri_tools.utils.files_are_unique( file_list = [ img_path.full( ".nii" )
	                                                         for img_path in img_paths
	                                                       ]
	                                         )
	      )

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

	# first, need to work out the slice acquisition times. These will be the same across all subjects,
	# but might as well extract one for each subject.
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


def mot_correct( conf, paths ):
	"""Performs motion correction"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running motion correction..." )

	st_paths = [ st_path.full( ".nii" ) for st_path in paths.func.sts ]
	corr_paths = [ corr_path.full( ".nii" ) for corr_path in paths.func.corrs ]

	# `mot_base` is one-based in the config
	i_base = conf.subj.mot_base - 1
	base_path = "{fname:s}[0]".format( fname = paths.func.sts[ i_base ].full( ".nii" ) )

	fmri_tools.preproc.mot_correct( orig_paths = st_paths,
	                                corr_paths = corr_paths,
	                                base_path = base_path,
	                                mc_path = paths.summ.motion.full( ".txt" ),
	                              )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( epi_paths = corr_paths,
	                                      out_path = paths.summ.corr.full()
	                                    )


def fieldmaps( conf, paths ):
	"""Prepare the fieldmaps"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running fieldmap preparation..." )

	# first, need to resample the mag and phase images to match the functionals

	# use the first corrected image as the image spec
	master = paths.func.corrs[ 0 ].full( ".nii[0]" )

	ref = conf.subj.subj_id + "_" + conf.exp.id + "-ref.nii"

	os.chdir( paths.fmap.base.full() )

	xtr_cmd = [ "3dcalc",
	            "-a", master,
	            "-expr", "a",
	            "-prefix", ref
	          ]

	fmri_tools.utils.run_cmd( " ".join( xtr_cmd ) )

	# use nipy to do the resampling; AFNI's doesn't seem to work properly
	ref_img = nipy.io.api.load_image( ref )

	for img_path in [ paths.fmap.mag.file( ".nii" ), paths.fmap.ph.file( ".nii" ) ]:

		source_img = nipy.io.api.load_image( img_path )

		if "ph" in img_path:

			r = np.arange( 1, 4097 )

			h = np.histogram( source_img._data.flatten(), r )

			diff = r[ np.argmax( h[ 0 ] ) ] - 2048

			source_img._data = np.mod( source_img._data - diff, 4096 )

			nipy.io.api.save_image( source_img, img_path )

			source_img = nipy.io.api.load_image( img_path )

		# use NN resampling to avoid out-of-range problems with interpolation
		resampled_img = nipy.algorithms.resample.resample_img2img( source_img,
		                                                           ref_img,
		                                                           order = 1
		                                                         )

		nipy.io.api.save_image( resampled_img, img_path )


	fmri_tools.preproc.make_fieldmap( mag_path = paths.fmap.mag.full(),
	                                  ph_path = paths.fmap.ph.full(),
	                                  fmap_path = paths.fmap.fmap.full(),
	                                  delta_te_ms = conf.acq.delta_te_ms
	                                )


def unwarp( conf, paths ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running distortion correction..." )

	for ( corr_path, uw_path ) in zip( paths.func.corrs, paths.func.uws ):

		fmri_tools.preproc.unwarp( epi_path = corr_path.full(),
		                           fmap_path = paths.fmap.fmap.full(),
		                           uw_path = uw_path.full(),
		                           dwell_ms = conf.acq.dwell_ms,
		                           uw_direction = conf.acq.ph_enc_dir,
		                           pass_nocheck = False
		                         )

	uw_files = [ uw.full() for uw in paths.func.uws ]

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( uw_files, paths.summ.uw.full() )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( uw_files, paths.summ.mean.full() )


def sess_reg( conf, paths ):
	"""Coregisters the session anatomical with its mean functional"""

	logger = logging.getLogger( __name__ )

	logger.info( "Running registration..." )

	if conf.subj.extra_al_params:
		extra_al_params = conf.subj.extra_al_params
	else:
		extra_al_params = None

	fmri_tools.preproc.img_reg( reg_dir = paths.reg.base.full(),
	                            base_file = paths.reg.mean.file( "+orig" ),
	                            mov_file = paths.reg.anat_ref.file( "+orig" ),
	                            extra_al_params = extra_al_params,
	                            epi_strip = "3dAutomask"
	                          )



def vol_to_surf( conf, paths ):
	"""Converts the functional volume-based images to SUMA surfaces."""

	logger = logging.getLogger( __name__ )
	logger.info( "Running volume to surface projection..." )

	start_dir = os.getcwd()

	for ( uw_file, surf_file, run_dir ) in zip( paths.func.uws, paths.func.surfs, paths.func.runs ):

		os.chdir( run_dir.full() )

		for hemi in [ "lh", "rh" ]:

			spec_file = paths.reg.spec.full( "_{hemi:s}.spec".format( hemi = hemi ) )
			spec_file = spec_file.replace( conf.subj.subj_id, conf.subj.subj_id.split( "_" )[ 0 ] )

			surf_cmd = [ "3dVol2Surf",
			             "-spec", spec_file,
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "15",
			             "-f_index", "nodes",
			             "-sv", paths.reg.anat_reg.full( "+orig" ),
			             "-grid_parent", uw_file.full( ".nii" ),
			             "-out_niml", surf_file.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( " ".join( surf_cmd ) )

	os.chdir( start_dir )


def design_prep( conf, paths ):
	"""Prepares the designs for GLM analysis"""

	# these are the text files we want to write that will contain the timing information for each
	# condition
	log_files = [ open( paths.ana.stim_times.full( "_{c:s}.txt".format( c = c ) ), "w" )
	              for c in conf.ana.conds
	            ]

	# loop over each experimental run
	for i_exp_run in xrange( conf.subj.n_exp_runs ):

		# check whether each condition has a valid event in this run - we need to mark it in its log
		# file if not
		cond_valid_evt = [ False ] * len( conf.ana.conds )

		# runtime log files stored as MATLAB mat files
		log_path = paths.logs.logs.full( "{r:d}.mat".format( r = i_exp_run + 1 ) )

		# load the log file
		log = scipy.io.loadmat( log_path, struct_as_record = False, squeeze_me = True )

		# get the onset time of the first trial, which is in arbitrary units
		# for some subjects, this is coded in a special field
		try:
			t0 = log[ "p" ].trial_onset[ 0 ]

		# but not all...
		except AttributeError:

			# in which case we need to rely on the first trial log
			t0 = log[ "data_mat" ][ 0, 3 ]

			# but we need to check that the first trial was recorded
			if log[ "data_mat" ][ 0, 1 ] != 1:
				raise ValueError( "No trial onset and no button press during first trial" )

		# data is keypresses x param, where param is:
		#  0 : subject number
		#  1 : trial number
		#  2 : response timestamp (same units as ``t0``)
		#  3 : stimulus onset (same units as ``t0``)
		#  4 : button one pressed (diamond)
		#  5 : button two pressed (non-diamond)
		#  6 : how long stimulus was up
		#  7 : how long its been since the last switch
		#  8 : trial type; 1 = short events (4s), 2 = longer (variable duration) events
		data = log[ "data_mat" ]

		# number of button presses; not necessarily same as the number of trials
		n_press = data.shape[ 0 ]

		# go through each button press
		for i_press in xrange( n_press ):

			# if we aren't on the last trial AND the next press has the same trial number, ignore this one
			# and go onto the next

			if ( ( i_press < ( n_press - 1 ) ) and
			     ( data[ i_press + 1, 1 ] == data[ i_press, 1 ] )
			   ):
				continue

			# subtract run start time offset to get into run-relative units
			trial_start = data[ i_press, 3 ] - t0

			if conf.subj.trig_on_ref:
				trial_start += conf.acq.tr_s

			# make sure they've pressed one, and only one, button
			assert( np.sum( data[ i_press, [ 4, 5 ] ] ) == 1 )

			sel_diamond = None

			# find which of the button indices (4 or 5) is True
			if data[ i_press, 4 ] == 1:
				sel_diamond = True
			elif data[ i_press, 5 ] == 1:
				sel_diamond = False

			if sel_diamond == None:
				raise ValueError( "Shouldnt happen" )

			# if a short event, duration is constant. Otherwise, it depends.
			if data[ i_press, 8 ] == 1:
				is_short = True
				duration = 4.0
			elif data[ i_press, 8 ] == 2:
				is_short = False
				duration = data[ i_press, 6 ]

			# get the onset in an AFNI friendly format
			time_str = "{o:.3f}".format( o = trial_start )

			# work out which condition is appropriate, based on trial type and button press
			if sel_diamond:
				if is_short:
					trial_cond = "di_short"
				else:
					# selected diamond to end the trial, meaning that most of it was non-diamond
					trial_cond = "nd_long"
			else:
				if is_short:
					trial_cond = "nd_short"
				else:
					# selected non-diamond to end the trial, meaning that most of it was diamond
					trial_cond = "di_long"

			# if its a 'long' trial, also need to add the duration info to the AFNI time string
			if trial_cond in [ "di_long", "nd_long" ]:
				time_str += ":{d:.3f}".format( d = duration )

			# find out the index corresponding to this condition type
			i_log = conf.ana.conds.index( trial_cond )

			# and write the trial info
			log_files[ i_log ].write( time_str + " " )

			# note that we've found a valid event for this condition
			cond_valid_evt[ i_log ] = True

		# end of this run, so now we need to see whether there were any missing conditions
		for ( i_log, cond_valid ) in enumerate( cond_valid_evt ):

			# mark any empty conditions by a '*' on that line in the file
			if not cond_valid:
				log_files[ i_log ].write( "*" )

		# end of the run, so write a newline to each log file
		_ = [ log_file.write( "\n" ) for log_file in log_files ]

	# all done
	_ = [ log_file.close() for log_file in log_files ]
