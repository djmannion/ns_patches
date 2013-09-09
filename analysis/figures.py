
import os, os.path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

import fmri_tools.utils

import ns_patches.config, ns_patches.paths


def _set_defaults():
	"""Set some sane defaults for figures.
	"""

	params = { 'axes.labelsize': 9 * ( 1 / 1.25 ),
	           'axes.titlesize' : 10,
	           'font.family' : 'Arial',
	           'font.sans-serif' : 'Helvetica',
	           'text.fontsize': 12,
	           'legend.fontsize': 7,
	           'xtick.labelsize': 8 * ( 1 / 1.25 ),
	           'xtick.direction' : 'out',
	           'xtick.major.size' : 2,
	           'ytick.labelsize': 8 * ( 1 / 1.25 ),
	           'ytick.direction' : 'out',
	           'ytick.major.size' : 2
	         }
	
	plt.rcParams.update( params )

	plt.ioff()


def _cleanup_fig( ax ):
	"""Apply some standard commands to clean up the axes on figures.
	"""

	for loc, spine in ax.spines.iteritems():

		spine.set_linewidth( 0.5 )

		if loc in [ "left", "bottom" ]:
			spine.set_position( ( "outward", 5 ) )
		elif loc in [ "right", "top" ]:
			spine.set_color( "none" )
		else:
			raise ValueError( "Unknown spine location: %s" % loc )

	ax.xaxis.set_ticks_position( "bottom" )
	ax.yaxis.set_ticks_position( "left" )


def all_subj_scatter():

	all_conf = ns_patches.config.get_conf( None, True )

	subj_ids = all_conf.all_subj.subj.keys()

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 7.08661, 10, forward = True )

	gs = gridspec.GridSpec( 3, 2 )

	for ( i_subj, subj_id ) in enumerate( subj_ids ):

		conf = ns_patches.config.get_conf( subj_id )

		paths = ns_patches.paths.get_subj_paths( conf )

		ax = plt.subplot( gs[ i_subj ] )

		ax.hold( True )

		# patches x ( coh, non-coh )
		data = np.loadtxt( paths.ana.vec_resp.full( ".txt" ) )

		# slope, intercept
		coef = np.loadtxt( paths.ana.regress.full( ".txt" ) )

		ax.scatter( data[ :, 0 ],
		            data[ :, 1 ],
		            facecolors = "None",
		            edgecolors = "k"
		          )

		_cleanup_fig( ax )

		xlim = plt.xlim()
		ylim = plt.ylim()

		max = np.max( [ xlim, ylim ] )
		min = np.min( [ xlim, ylim ] )

		# axis line
		ax.plot( [ min, max ], [ 0, 0 ], color = [ 0.5 ] * 3 )
		ax.plot( [ 0, 0 ], [ min, max ], color = [ 0.5 ] * 3 )

		# unity line
		ax.plot( [ min, max ], [ min, max ], "b--" )

		fit_x = np.linspace( min, max, 100 )

		lin_fit = np.polyval( coef, fit_x )

		ax.plot( fit_x, lin_fit, "g" )

		ax.set_xlim( [ min, max ] )
		ax.set_ylim( [ min, max ] )

		ax.text( 0.4, 0.05,
		         ", ".join( [ subj_id,
		                      "Slope: {s:.2f}".format( s = coef[ 0 ] ),
		                      "Intercept: {s:.2f}".format( s = coef[ 1 ] )
		                    ]
		                  ),
		         transform = ax.transAxes,
		         fontsize = 8 / 1.25
		       )

		if i_subj == len( subj_ids ) - 1:
			ax.set_xlabel( "Coherent response (psc)" )
			ax.set_ylabel( "Non-coherent response (psc)" )

	plt.subplots_adjust( left = 0.08,
	                     bottom = 0.06,
	                     right = 0.95,
	                     top = 0.97,
	                     wspace = 0.2,
	                     hspace = 0.2
	                   )

	plt.show()

