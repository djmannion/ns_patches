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

	os.chdir( paths.loc.base.full() )

	# minus one because the range is inclusive
	censor_vols = conf.loc.n_cull_vol - 1
	# in AFNI-aware format; (runs):start-end
	censor_str = "*:0-{v:.0f}".format( v = censor_vols )

	# model as a typical SPM event
	model_str = "'SPMG1({d:.0f})'".format( d = conf.loc.dur_s )

	for hemi in [ "lh", "rh" ]:

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ paths.func.surfs[ n - 1 ].full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		                  for n in conf.subj.loc_runs
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf.acq.tr_s ),
		                  "-polort", "a",  # auto baseline degree
		                  "-local_times",
		                  "-CENSORTR", censor_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "{n:d}".format( n = conf.stim.n_patches )
		                ]
		              )

		for i_patch in xrange( conf.stim.n_patches ):

			glm_cmd.extend( [ "-stim_times",
			                  "{n:d}".format( n = i_patch + 1 ),
			                  paths.loc.timing_base.full( "_{n:02d}.txt".format( n = i_patch ) ),
			                  model_str
			                ]
			              )

			glm_cmd.extend( [ "-stim_label",
			                  "{n:d}".format( n = i_patch + 1 ),
			                  "p{n:02d}".format( n = i_patch )
			                ]
			              )

		# run this first GLM
		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-Rbeta", paths.loc.beta.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		             "-tout",
		             "-Rbuck", paths.loc.glm.file( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( "'" +
		                 " ".join( [ paths.func.surfs[ n - 1 ].full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		                             for n in conf.subj.loc_runs
		                           ]
		                         ) +
		                 "'"
		               )

		# run the proper GLM
		fmri_tools.utils.run_cmd( " ".join( reml_cmd ) )

	os.chdir( start_dir )


def loc_patch_id( conf, paths ):
	"""Form a mask from the GLM output"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running localiser patch identification..." )

	# these are the t-statistics for each patch
	beta_bricks = "[2..64(2)]"

	# 3.3 == p of 0.001
	t_cutoff = 3.3 #2.58

	os.chdir( paths.loc.base.full() )

	for hemi in [ "lh", "rh" ]:

		glm_path = ( paths.loc.glm.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ) +
		             beta_bricks
		           )

		amax_cmd = [ "3dTstat",
		             "-overwrite",
		             "-argmax1",
		             "-max",
		             "-prefix", paths.loc.patch_id.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		             glm_path
		           ]

		fmri_tools.utils.run_cmd( " ".join( amax_cmd ) )

		thr_cmd = [ "3dcalc",
		            "-a", paths.loc.patch_id.full( "_{hemi:s}.niml.dset[1]".format( hemi = hemi ) ),
		            "-expr", "'step(a-{t:.4f})'".format( t = t_cutoff ),
		            "-prefix", paths.loc.patch_thr.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( " ".join( thr_cmd ) )

		bk_dset = paths.loc.patch.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )

		bk_cmd = [ "3dbucket",
		           "-overwrite",
		           "-prefix", bk_dset,
		           paths.loc.patch_id.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) ),
		           paths.loc.patch_thr.full( "_{hemi:s}.niml.dset".format( hemi = hemi ) )
		         ]

		fmri_tools.utils.run_cmd( " ".join( bk_cmd ) )

		full_bk_dset = paths.loc.patch.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		fmri_tools.utils.sparse_to_full( in_dset = bk_dset,
		                                 out_dset = full_bk_dset,
		                                 pad_node = "{n:d}".format( n = conf.subj.node_k[ hemi ] )
		                               )

		roi_path = paths.roi.vl.full( "_{hemi:s}-full.niml.dset".format( hemi = hemi ) )

		txt_path = paths.loc.patch.full( "_{hemi:s}.txt".format( hemi = hemi ) )

		# 3dmaskdump won't overwrite, so need to manually remove any previous file
		if os.path.exists( txt_path ):
			os.remove( txt_path )

		xtr_cmd = [ "3dmaskdump",
		            "-mask", full_bk_dset + "[2]",
		            "-noijk",
		            "-o", txt_path,
		            roi_path,
		            full_bk_dset + "[0]"
		          ]

		fmri_tools.utils.run_cmd( " ".join( xtr_cmd ) )



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
