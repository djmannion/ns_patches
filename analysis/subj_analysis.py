"""
Set of routines to analyse single-subject fMRI data for the diamond
 fMRI experiment.
"""

from __future__ import division

import os, os.path
import logging

import numpy as np

import fmri_tools.utils


def glm( conf, paths ):
	"""Experiment GLM"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running GLM..." )

	start_dir = os.getcwd()

	os.chdir( paths.ana.base.full() )

	n_cond = len( conf.ana.conds )

	cond_labels = conf.ana.conds

	cond_times = [ paths.ana.stim_times.full( "_{c:s}.txt".format( c = c ) )
	               for c in cond_labels
	             ]

	# minus one because the range is inclusive
	censor_vols = conf.ana.n_censor_vols - 1
	# in AFNI-aware format; (runs):start-end
	censor_str = "*:0-{v:.0f}".format( v = censor_vols )

	for hemi in [ "lh", "rh" ]:

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		surf_files = []

		for exp_run in conf.subj.exp_runs:
			surf_file = paths.func.surfs[ exp_run - 1 ]
			surf_files.append( surf_file.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ) )

		glm_cmd.extend( surf_files )

		glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf.acq.tr_s ),
		                  "-polort", "a",  # auto baseline degree
		                  "-local_times",
		                  "-CENSORTR", censor_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "{n:d}".format( n = n_cond )
		                ]
		              )

		# loop through each condition
		for i_cond in xrange( n_cond ):

			cond_label = cond_labels[ i_cond ]
			cond_time = cond_times[ i_cond ]
			cond_hrf = conf.ana.hrf_models[ i_cond ]
			cond_model = conf.ana.cond_models[ i_cond ]

			glm_cmd.extend( [ "-stim_label",
			                  "{sl:d}".format( sl = ( i_cond + 1 ) ),
			                  cond_label
			                ]
			              )

			glm_cmd.extend( [ "-" + cond_model,
			                  "{sl:d}".format( sl = ( i_cond + 1 ) ),
			                  cond_time,
			                  cond_hrf
			                ]
			              )

		# run this first GLM
		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-Rbeta", paths.ana.beta.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		             "-tout",
		             "-Rbuck", paths.ana.glm.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( "'" + " ".join( surf_files ) + "'" )

		# run the proper GLM
		fmri_tools.utils.run_cmd( " ".join( reml_cmd ) )

	os.chdir( start_dir )


def loc_glm( conf, paths ):
	"""Localiser GLMs"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running localiser GLM..." )

	start_dir = os.getcwd()

		loc_paths = paths.loc.loc[ loc_type ]

		os.chdir( loc_paths.base.full() )

		# minus one because the range is inclusive
		censor_vols = conf.ana.loc_n_censor_vols - 1
		# in AFNI-aware format; (runs):start-end
		censor_str = "*:0-{v:.0f}".format( v = censor_vols )

		# some subjects also have censoring at the end of the run
		if ( ( "stim_loc_end_censor" in conf.subj.__dict__ ) and
		     ( conf.subj.stim_loc_end_censor ) and
		     ( loc_type == "stim" )
		   ):

			# special case for this subject with more volumes than expected
			censor_str += " *:127..134"

		# need to calculate the rows in the motion correction file that correspond to this localiser runs
		run_num = conf.subj.loc_runs[ i_loc ]
		start_vol = np.sum( conf.subj.run_n_vols[ :( run_num - 1 ) ] )
		end_vol = start_vol + conf.subj.run_n_vols[ run_num - 1 ] - 1

		assert( np.arange( start_vol, end_vol + 1 ).shape[ 0 ] == conf.subj.run_n_vols[ run_num - 1 ] )

		# timing of the onset of 'ON' blocks
		block_starts = np.arange( conf.ana.loc_n_vol_per_block,
		                          conf.subj.run_n_vols[ run_num - 1 ],
		                          conf.ana.loc_n_vol_per_block * 2  # doubled because we only model half
		                        )
		# convert to seconds
		block_starts *= conf.acq.tr_s

		if conf.subj.trig_on_ref:
			block_starts += conf.acq.tr_s

		# and put into AFNI format
		time_str = "'1D: " + " ".join( [ "{t:.3f}".format( t = t ) for t in block_starts ] ) + "'"

		# model as a typical SPM block
		model_str = "SPMG1({d:.0f})".format( d = conf.ana.loc_n_vol_per_block * conf.acq.tr_s )

		for hemi in [ "lh", "rh" ]:

			glm_cmd = [ "3dDeconvolve",
			            "-input",
			            paths.func.surfs[ run_num - 1 ].full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
			          ]

			glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf.acq.tr_s ),
			                  "-polort", "a",  # auto baseline degree
			                  "-local_times",
			                  "-CENSORTR", censor_str,
			                  "-xjpeg", "exp_design.png",
			                  "-x1D", "exp_design",
			                  "-overwrite",
			                  "-x1D_stop",  # want to use REML, so don't bother running
			                  "-num_stimts", "1",
			                  "-stim_label", "1", "ON",
			                  "-stim_times", "1", time_str, model_str
			                ]
			              )

			# run this first GLM
			fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

			# delete the annoying command file that 3dDeconvolve writes
			os.remove( "Decon.REML_cmd" )

			reml_cmd = [ "3dREMLfit",
			             "-matrix", "exp_design.xmat.1D",
			             "-Rbeta", loc_paths.beta.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
			             "-tout",
			             "-Rbuck", loc_paths.glm.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
			             "-overwrite",
			             "-input",
			             paths.func.surfs[ run_num - 1 ].full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
			           ]

			# run the proper GLM
			fmri_tools.utils.run_cmd( " ".join( reml_cmd ) )

	os.chdir( start_dir )


def loc_mask( conf, paths ):
	"""Form a mask from the GLM output"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running localising mask creation..." )

	# these are the indices into the GLM files for the data we want to convert
	loc_brick = "[2]"

	loc_paths = paths.loc.loc[ "stim" ]

	for hemi in [ "lh", "rh" ]:

		glm_path = ( loc_paths.glm.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ) +
		             loc_brick
		           )

		# check that the label is as expected
		glm_label = fmri_tools.utils.get_dset_label( glm_path )
		assert( glm_label == [ "ON#0_Tstat" ] )

		fdr_path = loc_paths.fdr.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		mask_path = loc_paths.mask.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

		fmri_tools.utils.loc_mask( glm_path = glm_path,
		                           fdr_path = fdr_path,
		                           mask_path = mask_path,
		                           q_thresh = conf.ana.loc_q,
		                           pos_only = False
		                         )

		pad_mask_path = loc_paths.mask.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )
		pad_nodes = "{nk:d}".format( nk = conf.subj.node_k[ hemi ] )

		fmri_tools.utils.sparse_to_full( in_dset = mask_path,
		                                 out_dset = pad_mask_path,
		                                 pad_node = pad_nodes
		                               )


def beta_to_psc( conf, paths ):
	"""Convert the GLM beta weights into units of percent signal change"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running beta to PSC conversion..." )

	for hemi in [ "lh", "rh" ]:

		beta_path = paths.ana.beta.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

		# this returns the labels of all bricks in the dataset
		beta_labels = fmri_tools.utils.get_dset_label( beta_path )

		beta_bricks = []

		for ( i_beta, beta_label ) in enumerate( beta_labels ):

			if beta_label.split( "#" )[ 0 ] in conf.ana.conds:
				beta_bricks.append( "{n:d}".format( n = i_beta ) )

		beta_bricks = "[" + ",".join( beta_bricks ) + "]"

		design_path = ( paths.ana.base + "exp_design.xmat.1D" ).full()

		bltc_path = paths.ana.bltc.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		bl_path = paths.ana.bl.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		psc_path = paths.ana.psc.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

		fmri_tools.utils.beta_to_psc( beta_path = beta_path,
		                              beta_bricks = beta_bricks,
		                              design_path = design_path,
		                              bltc_path = bltc_path,
		                              bl_path = bl_path,
		                              psc_path = psc_path
		                            )

		pad_psc_path = paths.ana.psc.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )
		pad_nodes = "{nk:d}".format( nk = conf.subj.node_k[ hemi ] )

		fmri_tools.utils.sparse_to_full( in_dset = psc_path,
		                                 out_dset = pad_psc_path,
		                                 pad_node = pad_nodes
		                               )

	# now do the same thing for the localisers

	# these are the indices into the beta files for the data we want to convert
	beta_bricks = "[3]"

	for loc_paths in paths.loc.loc.values():

		for hemi in [ "lh", "rh" ]:

			beta_path = loc_paths.beta.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

			design_path = ( loc_paths.base + "exp_design.xmat.1D" ).full()

			bltc_path = loc_paths.bltc.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
			bl_path = loc_paths.bl.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
			psc_path = loc_paths.psc.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

			fmri_tools.utils.beta_to_psc( beta_path = beta_path,
			                              beta_bricks = beta_bricks,
			                              design_path = design_path,
			                              bltc_path = bltc_path,
			                              bl_path = bl_path,
			                              psc_path = psc_path
			                            )

			pad_psc_path = loc_paths.psc.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )
			pad_nodes = "{nk:d}".format( nk = conf.subj.node_k[ hemi ] )

			fmri_tools.utils.sparse_to_full( in_dset = psc_path,
			                                 out_dset = pad_psc_path,
			                                 pad_node = pad_nodes
			                               )


def roi_xtr( conf, paths ):
	"""Extract PSC and statistics data from ROIs"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running ROI extraction..." )

	start_dir = os.getcwd()

	os.chdir( paths.roi.base.dir() )

	for hemi in [ "lh", "rh" ]:

		roi_path = paths.roi.rois.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		dset_path = paths.ana.psc.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		loc_path = paths.loc.loc[ "stim" ].psc.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		mask_path = paths.loc.loc[ "stim" ].mask.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		txt_path = paths.roi.psc.full( "_{hemi:s}.txt".format( hemi = hemi ) )

		# 3dmaskdump won't overwrite, so need to manually remove any previous file
		if os.path.exists( txt_path ):
			os.remove( txt_path )

		xtr_cmd = [ "3dmaskdump",
		            "-mask", roi_path,
		            "-noijk",
		            "-o", txt_path,
		            roi_path,
		            dset_path,
		            loc_path,
		            mask_path
		          ]

		fmri_tools.utils.run_cmd( " ".join( xtr_cmd ) )

		# write out the column headers to a different file
		headers = []

		beta_labels = fmri_tools.utils.get_dset_label( dset_path )

		headers.append( "ROI" )
		headers.extend( beta_labels )
		headers.append( "Localiser PSC" )
		headers.append( "Localiser mask" )

		with open( paths.roi.psc_header.full( ".txt" ), "w" ) as header_file:
			header_file.write( "\t".join( headers ) )

	os.chdir( start_dir )
