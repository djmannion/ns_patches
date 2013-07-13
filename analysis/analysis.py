"""Performs the analysis for a single subject"""

import os, os.path
import logging

import numpy as np

import fmri_tools.utils


def _get_timing( conf, paths ):

	timing = []

	pre_len_s = conf.exp.n_censor_vols * conf.acq.tr_s

	for run_num in xrange( 1, conf.subj.n_runs + 1 ):

		# load the run log file
		run_ext = "{r:02d}_log.npz".format( r = run_num )
		run_log = np.load( paths.logs.run_log_base.full( run_ext ) )

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

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		surf_paths = [ surf_path.full( "_{h:s}.niml.dset".format( h = hemi ) )
		               for surf_path in paths.func.surfs
		             ]

		glm_cmd.extend( surf_paths )

		glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf.acq.tr_s ),
		                  "-polort", "a",  # auto baseline degree
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
			                  "'1D: " + reg_times + "'",
			                  model_str
			                ]
			              )

		# run this first GLM
#		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		print " ".join( glm_cmd )

		return

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


