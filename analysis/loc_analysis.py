"""Performs the localiser analysis for a single subject"""

import os, os.path
import logging

import scipy.stats

import fmri_tools.utils


def glm( conf, paths ):
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

		surf_paths = [ surf_path.full( "_{h:s}.niml.dset".format( h = hemi ) )
		               for surf_path in paths.func.surfs
		             ]

		glm_cmd.extend( surf_paths )

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

			timing_ext = "_{n:02d}.txt".format( n = i_patch )

			timing_file = paths.loc.timing_base.full( timing_ext )

			glm_cmd.extend( [ "-stim_times",
			                  "{n:d}".format( n = i_patch + 1 ),
			                  timing_file,
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

		beta_file = paths.loc.beta.file( "_{h:s}.niml.dset".format( h = hemi ) )
		buck_file = paths.loc.glm.file( "_{h:s}.niml.dset".format( h = hemi ) )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
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


def patch_id( conf, paths ):
	"""Form a mask from the GLM output"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running localiser patch identification..." )

	# these are the t-statistics for each patch
	t_bricks = "[2..64(2)]"

	t_cutoff = scipy.stats.t.isf( conf.loc.t_p, conf.loc.dof )

	os.chdir( paths.loc.base.full() )

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "_{h:s}.niml.dset".format( h = hemi )

		glm_path = paths.loc.glm.full( hemi_ext + t_bricks )

		# first, mark each regressor as significant or not
		sig_path = paths.loc.sig.full( hemi_ext )

		sig_cmd = [ "3dcalc",
		            "-a", glm_path,
		            "-expr", "ispositive(step(a-{t:.4f}))".format( t = t_cutoff ),
		            "-prefix", sig_path,
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( " ".join( sig_cmd ) )


		# then, count how many significant regressors there are at each node
		sig_sum_path = paths.log.sig_sum.full( hemi_ext )

		sig_sum_cmd = [ "3dTstat",
		                "-overwrite",
		                "-sum",
		                "-prefix", sig_sum_path,
		                sig_path
		              ]

		fmri_tools.utils.run_cmd( " ".join( sig_sum_cmd ) )


		# now work out which ID is significant for each node, subject to the
		# constraint that there is only one significant patch
		id_path = paths.loc.patch_id.full( hemi_ext )

		id_cmd = [ "3dTstat",
		           "-overwrite",
		           "-argmax1",
		           "-prefix", id_path,
		           "-mask", sig_sum_path,
		           "-mrange", "1", "1",
		           glm_path
		         ]

		fmri_tools.utils.run_cmd( " ".join( id_cmd ) )

		full_hemi_ext = "_{h:s}-full.niml.dset".format( h = hemi )

		# need to pad to full for integration with ROIs
		id_path_full = paths.loc.patch_id.full( full_hemi_ext )

		pad_k = "{n:d}".format( n = conf.subj.node_k[ hemi ] )

		fmri_tools.utils.sparse_to_full( in_dset = id_path,
		                                 out_dset = id_path_full,
		                                 pad_node = pad_k
		                               )

		roi_path = paths.roi.vl.full( full_hemi_ext )

		txt_path = paths.loc.patch.full( hemi_ext )

		# 3dmaskdump won't overwrite, so need to manually remove any previous file
		if os.path.exists( txt_path ):
			os.remove( txt_path )

		xtr_cmd = [ "3dmaskdump",
		            "-noijk",
		            "-o", txt_path,
		            roi_path,
		            id_path_full
		          ]

		fmri_tools.utils.run_cmd( " ".join( xtr_cmd ) )
