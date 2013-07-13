"""Performs the analysis for a single subject"""

import os, os.path
import logging

import numpy as np

import fmri_tools.utils


def glm( conf, paths ):
	"""Experiment GLM"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running GLM..." )

	start_dir = os.getcwd()

	# write timing files
	timing_path = paths.ana.stim_times.full( ".txt" )

	with open( timing_path, "w" ) as timing_file:

		for run_num in xrange( 1, conf.subj.n_runs + 1 ):

			# load the run log file
			run_ext = "_{r:d}_log.npz".format( r = run_num )
			run_log = np.load( paths.logs.run_log_base.full( run_ext ) )

			# we're interested in the 'seq' variable
			run_seq = run_log[ "seq" ]

			pre_len_s = conf.exp.n_censor_vols * conf.acq.tr_s

			# estimate 26s as the duration of an event
			# if we don't cull these events, they are zero in the GLM and that messes
			# stuff up
			cull_from_pre_cutoff = pre_len_s - 26

			i_keep = ( run_seq > cull_from_pre_cutoff )

			run_seq = run_seq[ i_keep ]

			# convert to strings
			seq_str = [ "{n:.02f}".format( n = n ) for n in run_seq ]

			timing_file.write( " ".join( seq_str ) + "\n" )

	os.chdir( paths.ana.base.full() )

	# minus one because the range is inclusive
	censor_vols = conf.exp.n_censor_vols - 1
	# in AFNI-aware format; (runs):start-end
	censor_str = "*:0-{v:.0f}".format( v = censor_vols )

	# model as a typical SPM event
	model_str = "'SPMG1({d:.0f})'".format( d = conf.exp.img_on_s )

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
		                  "-num_stimts", "1",
		                  "-stim_times_IM", "1", timing_path, model_str,
		                  "-stim_label", "1", "stim"
		                ]
		              )

		# run this first GLM
		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		beta_file = paths.ana.beta.file( "_{h:s}.niml.dset".format( h = hemi ) )
		buck_file = paths.ana.glm.file( "_{h:s}.niml.dset".format( h = hemi ) )

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


