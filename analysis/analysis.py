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

		con_str = ( "SYM: " +
		            " ".join( [ "+r{x:s}".format( x = x )
		                        for x in regressors[ 0 ]
		                      ]
		                    )
		          )

		glm_cmd.extend( [ "-gltsym", "'" + con_str + "'",
		                  "-glt_label", "1", "all"
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

	control_resp = np.empty( img_resp.shape )
	control_resp.fill( np.NAN )

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

		# find the coherent image for each trial
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

				np.random.shuffle( i_coh_trials )
				np.random.shuffle( i_incoh_trials )

				p_one_mean = np.mean( [ patch_data[ i_patch, i_coh_trials[ 0 ] ],
				                        patch_data[ i_patch, i_incoh_trials[ 0 ] ]
				                      ]
				                    )
				control_resp[ i_patch, i_img, run_num - 1, 0 ] = p_one_mean

				p_two_mean = np.mean( [ patch_data[ i_patch, i_coh_trials[ 1 ] ],
				                        patch_data[ i_patch, i_incoh_trials[ 1 ] ]
				                      ]
				                    )
				control_resp[ i_patch, i_img, run_num - 1, 1 ] = p_two_mean

	assert np.sum( np.isnan( img_resp ) ) == 0
	assert np.sum( np.isnan( control_resp ) ) == 0

	np.save( paths.ana.img_resp.full( ".npy" ), img_resp )
	np.save( paths.ana.img_resp.full( "-control.npy" ), control_resp )


def vec_resp( conf, paths ):

	# is 26 x 20 x 8 x 2
	img_resp = np.load( paths.ana.img_resp.full( ".npy" ) )

	# average over runs
	img_resp = np.mean( img_resp, axis = 2 )

	# concatenate patches and images
	img_resp = np.vstack( [ img_resp[ ..., 0 ].flatten(),
	                        img_resp[ ..., 1 ].flatten()
	                      ]
	                    )

	assert img_resp.shape[ 0 ] == 2

	img_resp = img_resp.T

	np.savetxt( paths.ana.vec_resp.full( ".txt" ), img_resp )


def subj_regress( conf, paths ):

	vec_resp = np.loadtxt( paths.ana.vec_resp.full( ".txt" ) )

	reg_coef = scipy.stats.linregress( vec_resp[ :, 0 ],
	                                   vec_resp[ :, 1 ]
	                                 )[ :2 ]

	reg_coef = total_least_squares( vec_resp[ :, 0 ],
	                                vec_resp[ :, 1 ]
	                              )

	np.savetxt( paths.ana.regress.full( ".txt" ), reg_coef )


def total_least_squares(data1, data2, data1err=None, data2err=None,
        print_results=False, ignore_nans=True, intercept=True,
        return_error=False, inf=1e10):
    """
Use Singular Value Decomposition to determine the Total Least Squares linear fit to the data.
(e.g. http://en.wikipedia.org/wiki/Total_least_squares)
data1 - x array
data2 - y array

if intercept:
returns m,b in the equation y = m x + b
else:
returns m

print tells you some information about what fraction of the variance is accounted for

ignore_nans will remove NAN values from BOTH arrays before computing

Parameters
----------
data1,data2 : np.ndarray
Vectors of the same length indicating the 'x' and 'y' vectors to fit
data1err,data2err : np.ndarray or None
Vectors of the same length as data1,data2 holding the 1-sigma error values

"""

    if ignore_nans:
        badvals = np.isnan(data1) + np.isnan(data2)
        if data1err is not None:
            badvals += np.isnan(data1err)
        if data2err is not None:
            badvals += np.isnan(data2err)
        goodvals = True-badvals
        if goodvals.sum() < 2:
            if intercept:
                return 0,0
            else:
                return 0
        if badvals.sum():
            data1 = data1[goodvals]
            data2 = data2[goodvals]

    
    if intercept:
        dm1 = data1.mean()
        dm2 = data2.mean()
    else:
        dm1,dm2 = 0,0

    arr = np.array([data1-dm1,data2-dm2]).T

    U,S,V = np.linalg.svd(arr, full_matrices=False)

    # v should be sorted.
    # this solution should be equivalent to v[1,0] / -v[1,1]
    # but I'm using this: http://stackoverflow.com/questions/5879986/pseudo-inverse-of-sparse-matrix-in-python
    M = V[-1,0]/-V[-1,-1]

    varfrac = S[0]/S.sum()*100
    if varfrac < 50:
        raise ValueError("ERROR: SVD/TLS Linear Fit accounts for less than half the variance; this is impossible by definition.")

    # this is performed after so that TLS gives a "guess"
    if data1err is not None or data2err is not None:
        try:
            from scipy.odr import RealData,Model,ODR
        except ImportError:
            raise ImportError("Could not import scipy; cannot run Total Least Squares")

        def linmodel(B,x):
            if intercept:
                return B[0]*x + B[1]
            else:
                return B[0]*x

        if data1err is not None:
            data1err = data1err[goodvals]
            data1err[data1err<=0] = inf
        if data2err is not None:
            data2err = data2err[goodvals]
            data2err[data2err<=0] = inf

        if any([data1.shape != other.shape for other in (data2,data1err,data2err)]):
            raise ValueError("Data shapes do not match")

        linear = Model(linmodel)
        data = RealData(data1,data2,sx=data1err,sy=data2err)
        B = data2.mean() - M*data1.mean()
        beta0 = [M,B] if intercept else [M]
        myodr = ODR(data,linear,beta0=beta0)
        output = myodr.run()

        if print_results:
            output.pprint()

        if return_error:
            return np.concatenate([output.beta,output.sd_beta])
        else:
            return output.beta



    if intercept:
        B = data2.mean() - M*data1.mean()
        if print_results:
            print "TLS Best fit y = %g x + %g" % (M,B)
            print "The fit accounts for %0.3g%% of the variance." % (varfrac)
            print "Chi^2 = %g, N = %i" % (((data2-(data1*M+B))**2).sum(),data1.shape[0]-2)
        return M,B
    else:
        if print_results:
            print "TLS Best fit y = %g x" % (M)
            print "The fit accounts for %0.3g%% of the variance." % (varfrac)
            print "Chi^2 = %g, N = %i" % (((data2-(data1*M))**2).sum(),data1.shape[0]-1)
        return M




