"""Performs the localiser analysis for a single subject"""

import os, os.path
import logging

import fmri_tools.utils


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

		glm_path = ( paths.loc.glm.full( "_{h:s}.niml.dset".format( h = hemi ) ) +
		             beta_bricks
		           )

		amx_file = paths.loc.patch_id.full( "_{h:s}.niml.dset".format( h = hemi ) )

		amax_cmd = [ "3dTstat",
		             "-overwrite",
		             "-argmax1",
		             "-max",
		             "-prefix", amx_file,
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
