"""Performs the analysis for a single subject"""

import os, os.path
import logging

import numpy as np
import scipy.stats

import fmri_tools.utils

import ns_patches.config, ns_patches.paths


def _get_timing( conf, paths ):
	"Calculate the timing of each regressor in each run"

	timing = []

	pre_len_s = conf.exp.n_censor_vols * conf.acq.tr_s

	for run_num in xrange( 1, conf.subj.n_runs + 1 ):

		# load the run log file
		run_log = _load_run_log( run_num, paths )

		run_seq = run_log[ "seq" ]

		run_evt = []

		for evt_t in run_seq:

			# postpone pre events, for now
			if evt_t < pre_len_s:
				continue

			run_evt.append( [ evt_t ] )

		assert len( run_evt ) == conf.exp.n_trials

		# now go back to the pre events
		i_pre_evts = np.where( run_seq < pre_len_s )[ 0 ]

		for ( i_countback, i_pre ) in enumerate( i_pre_evts[ ::-1 ] ):
			run_evt[ -( i_countback + 1 ) ].append( run_seq[ i_pre ] )

		timing.append( run_evt )

	regress_num = 1

	regress = []

	for i_run in xrange( len( timing ) ):

		run_timing = timing[ i_run ]

		for i_evt in xrange( len( run_timing ) ):

			evt_global_s = [ t + ( i_run * conf.exp.run_len_s )
			                 for t in sorted( run_timing[ i_evt ] )
			               ]

			evt_str = [ "{n:.0f}".format( n = n ) for n in evt_global_s ]

			regress.append( [ "{n:d}".format( n = regress_num ),
			                  " ".join( evt_str )
			                ]
			              )

			regress_num += 1

	return regress


def glm( conf, paths ):
	"""Experiment GLM"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running GLM..." )

	start_dir = os.getcwd()

	regressors = _get_timing( conf, paths )

	os.chdir( paths.ana.base.full() )

	# minus one because the range is inclusive
	censor_vols = conf.exp.n_censor_vols - 1
	# in AFNI-aware format; (runs):start-end
	censor_str = "*:0-{v:.0f}".format( v = censor_vols )

	# model as a typical SPM event
	model_str = "'SPMG1({d:.0f})'".format( d = conf.exp.img_on_s )

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "_{h:s}".format( h = hemi )

		# write the ROI mask
		mask_file = paths.ana.mask.full( hemi_ext + "-full.niml.dset" )
		mask_cmd = [ "3dcalc",
		             "-a",
		             paths.ana.vl.full( hemi_ext + "-full.niml.dset" ),
		             "-expr", "equals(a,1)",
		             "-prefix", mask_file,
		             "-overwrite"
		           ]

		fmri_tools.utils.run_cmd( " ".join( mask_cmd ) )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		surf_paths = [ surf_path.full( "_{h:s}-full.niml.dset".format( h = hemi ) )
		               for surf_path in paths.func.surfs
		             ]

		glm_cmd.extend( surf_paths )

		glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf.acq.tr_s ),
		                  "-polort", "a",  # auto baseline degree
		                  "-mask", mask_file,
		                  "-global_times",
		                  "-CENSORTR", censor_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "{n:d}".format( n = len( regressors ) )
		                ]
		              )

		# now comes the behemoth
		for ( reg_num, reg_times ) in regressors:

			glm_cmd.extend( [ "-stim_label", reg_num, "r" + reg_num ] )

			glm_cmd.extend( [ "-stim_times",
			                  reg_num,
			                  "'1D: " + reg_times + "'",
			                  model_str
			                ]
			              )

		# run this first GLM
		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		beta_file = paths.ana.beta.file( hemi_ext + "-full.niml.dset" )
		buck_file = paths.ana.glm.file( hemi_ext + "-full.niml.dset" )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-mask", mask_file,
		             "-Rbeta", beta_file,
		             "-tout",
		             "-Rbuck", buck_file,
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( "'" + " ".join( surf_paths ) + "'" )

		# run the proper GLM
		fmri_tools.utils.run_cmd( " ".join( reml_cmd ) )

	os.chdir( start_dir )


def beta_to_psc( conf, paths ):
	"""Convert the beta estimates into units of percent signal change"""

	start_dir = os.getcwd()
	os.chdir( paths.ana.base.full() )

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "_{h:s}-full.niml.dset".format( h = hemi )

		beta_path = paths.ana.beta.full( hemi_ext )
		design_path = "exp_design.xmat.1D"

		# to write
		bltc_path = paths.ana.bltc.file( hemi_ext )
		bl_path = paths.ana.bl.file( hemi_ext )
		psc_path = paths.ana.psc.file( hemi_ext )

		# 4 orthogonal polynomial regressors per run
		n_nuisance = conf.subj.n_runs * 4


		beta_bricks = "[{n:d}..$]".format( n = n_nuisance )

		fmri_tools.utils.beta_to_psc( beta_path,
		                              beta_bricks,
		                              design_path,
		                              bltc_path,
		                              bl_path,
		                              psc_path,
		                            )

	os.chdir( start_dir )


def patch_dump( conf, paths ):
	"""Write out the PSC and patch ID for each node"""

	loc_id = conf.subj.subj_id + "_loc"
	loc_conf = ns_patches.config.get_conf( loc_id )
	loc_paths = ns_patches.paths.get_subj_paths( loc_conf )

	resp = []

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "_{h:s}".format( h = hemi )

		# location of the patch ID
		id_path = loc_paths.loc.patch_id.full( hemi_ext + "-full.niml.dset" )

		psc_path = paths.ana.psc.full( hemi_ext + "-full.niml.dset" )

		resp_path = paths.ana.patch_resp.full( hemi_ext + ".txt" )

		if os.path.exists( resp_path ):
			os.remove( resp_path )

		cmd = [ "3dmaskdump",
		        "-noijk",
		        "-mask", id_path,
		        "-o", resp_path,
		        id_path,
		        psc_path
		      ]

		fmri_tools.utils.run_cmd( " ".join( cmd ) )

		# we want to combine across hemispheres, so load the result
		resp.append( np.loadtxt( resp_path ) )

	# combine across hemispheres
	resp = np.vstack( resp )

	resp_path = paths.ana.patch_resp.full( ".txt" )

	np.savetxt( resp_path, resp )


def _load_run_log( run_num, paths ):
	"Loads the log file for a particular subject and run"

	# load the run log file
	run_ext = "{r:02d}_log.npz".format( r = run_num )
	run_log = np.load( paths.logs.run_log_base.full( run_ext ) )

	return run_log


def image_resp( conf, paths ):
	"""Calculate the response to each image"""

	resp_path = paths.ana.patch_resp.full( ".txt" )

	# nodes x trials
	resp = np.loadtxt( resp_path )

	# first step is to convert to patch x trial
	# the '-1' is because `resp` has a column for the patch id
	patch_data = np.empty( ( conf.exp.n_mod_patches, resp.shape[ 1 ] - 1 ) )
	patch_data.fill( np.NAN )

	for i_patch in xrange( conf.exp.n_mod_patches ):

		# in `resp`, patch ID is stored as one-based
		i_node_patch = np.where( resp[ :, 0 ] == ( i_patch + 1 ) )[ 0 ]

		# average over nodes for this patch
		patch_data[ i_patch, : ] = np.mean( resp[ i_node_patch, 1: ], axis = 0 )

	assert np.sum( np.isnan( patch_data ) ) == 0

	img_resp = np.empty( ( conf.exp.n_mod_patches,
	                       conf.exp.n_img,
	                       conf.subj.n_runs,
	                       2
	                     )
	                   )
	img_resp.fill( np.NAN )

	# now we need to load the run info
	for run_num in xrange( 1, conf.subj.n_runs + 1 ):

		run_log = _load_run_log( run_num, paths )

		# `img_trials` is patches x trials
		img_trials = run_log[ "img_trials" ]

		# cull the pre-trials
		img_trials = img_trials[ :, -conf.exp.n_trials: ]

		assert img_trials.shape[ 1 ] == conf.exp.n_trials

		# and limit to the modulated patches
		img_trials = img_trials[ conf.exp.mod_patches, : ]

		assert img_trials.shape[ 0 ] == conf.exp.n_mod_patches

		i_coh_img = scipy.stats.mode( img_trials, axis = 0 )[ 0 ][ 0 ]

		i_run_offset = ( run_num - 1 ) * conf.exp.n_trials

		for i_patch in xrange( conf.exp.n_mod_patches ):

			patch_info = img_trials[ i_patch, : ]

			for i_img in xrange( conf.exp.n_img ):

				# find out the trials where this patch showed the image
				i_coh_trials = np.where( np.logical_and( patch_info == i_img,
				                                         i_img == i_coh_img
				                                       )
				                       )[ 0 ]

				i_incoh_trials = np.where( np.logical_and( patch_info == i_img,
				                                           i_img != i_coh_img
				                                         )
				                         )[ 0 ]

				assert len( np.intersect1d( i_coh_trials, i_incoh_trials ) ) == 0

				assert len( i_coh_trials ) == 2
				assert len( i_incoh_trials ) == 2

				i_coh_trials += i_run_offset
				i_incoh_trials += i_run_offset

				coh_mean = np.mean( patch_data[ i_patch, i_coh_trials ] )
				img_resp[ i_patch, i_img, run_num - 1, 0 ] = coh_mean

				incoh_mean = np.mean( patch_data[ i_patch, i_incoh_trials ] )
				img_resp[ i_patch, i_img, run_num - 1, 1 ] = incoh_mean

	assert np.sum( np.isnan( img_resp ) ) == 0

	return img_resp
