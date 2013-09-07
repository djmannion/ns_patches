
import numpy as np
import scipy.stats

import ns_patches.config, ns_patches.paths


def reg_coeff_t( conf, paths ):

	subj_ids = conf.all_subj.subj.keys()

	n_subj = len( subj_ids )

	reg_coef = np.empty( ( n_subj, 2 ) )
	reg_coef.fill( np.NAN )

	for ( i_subj, subj_id ) in enumerate( subj_ids ):

		subj_conf = ns_patches.config.get_conf( subj_id )

		subj_paths = ns_patches.paths.get_subj_paths( subj_conf )

		subj_coef = np.loadtxt( subj_paths.ana.regress.full( ".txt" ) )

		reg_coef[ i_subj, : ] = subj_coef

	# slope t-test
	( slope_t, slope_p ) = scipy.stats.ttest_1samp( a = reg_coeff[ :, 0 ],
	                                                popmean = 1.0
	                                              )

	# intercept t-test
	( intrcpt_t, intrcpt_p ) = scipy.stats.ttest_1samp( a = reg_coef[ :, 1 ],
	                                                    popmean = 0.0 )

	return reg_coef
