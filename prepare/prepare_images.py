
import os, os.path
import numpy as np

import scipy.stats

import psychopy.filters

import stimuli.utils, stimuli.upenn_db

import ns_patches.config, ns_patches.paths, ns_patches.stimulus
import ns_patches.exp.exp


def img_browser( skip_discarded = True,
                 save_changes = False,
                 seek_to_first_A = False,
                 seek_to_id = None,
                 mon_name = "UMN_7T_colour",
                 shuffle = False
               ):
    """Interface to browse through and select from the image database.

    Parameters
    ----------
    skip_discarded : bool, optional
        Whether to show images marked as discarded.
    save_changes : bool, optional
        Whether to save any changes to the master file.
    seek_to_first_A : bool, optional
        Start at the first image marked 'A'.
    seek_to_id : int or None, optional
        Start at a particular ID

    """

    conf = ns_patches.config.get_conf()
    paths = ns_patches.paths.get_exp_paths( conf )

    img_db_info = read_img_db_info( paths.img_db_info )

    if shuffle:
        img_db_info = img_db_info[ np.random.permutation( len( img_db_info ) ) ]

    win = psychopy.visual.Window( ( 1024, 768 ),
                                  fullscr = False,
                                  allowGUI = True,
                                  monitor = mon_name,
                                  units = "pix"
                                )

    fix_stim = ns_patches.stimulus.get_fixation( win, conf )

    status = psychopy.visual.TextStim( win = win,
                                       text = "placeholder",
                                       pos = ( 0, -300 ),
                                       units = "pix",
                                       alignHoriz = "center"
                                     )

    i_first = 0

    try:
        if seek_to_first_A:
            i_first = np.where( img_db_info[ : ][ "status" ] == "A" )[ 0 ][ 0 ]

    except IndexError:
        win.close()
        print "No unclassified image found"
        return

    if seek_to_id is not None:
        i_first = np.where( img_db_info[ "id" ] == seek_to_id )[ 0 ][ 0 ]

    i_img = i_first

    # calculate a set of incoherent patches for each image
#   incoh = ns_patches.config.sel_incoh( n_patch = conf.exp.n_mod_patches,
#                                        n_incoh = conf.exp.n_incoh_patches,
#                                        n_trials = 80
#                                      )

    img_trials = ns_patches.config.gen_exp_patch_timing( conf )

    incoh = np.zeros( ( conf.exp.n_incoh_patches, conf.exp.n_trials ) )

    for i_trial in xrange( incoh.shape[ 1 ] ):
        i_coh = scipy.stats.mode( img_trials[ i_trial, : ] )[ 0 ]

        i_incoh = np.where( img_trials[ i_trial, : ] != i_coh )[ 0 ]

        incoh[ :, i_trial ] = i_incoh

    masks = ns_patches.exp.exp.gen_masks( conf )
    stim = ns_patches.exp.exp.gen_stim( conf, win, masks )

    change_image = False

    i_trial = 0

    keep_going = True

    while keep_going:

        incoh_patches = [ conf.exp.mod_patches[ int( i ) ] for i in incoh[ :, i_trial ] ]
        coh_patches = np.setdiff1d( np.arange( conf.stim.n_patches ), incoh_patches )

        img_ok = False

        while not img_ok:

            if i_img == len( img_db_info ):
                i_img = 0

            coh_img_info = img_db_info[ i_img ]

            if skip_discarded and coh_img_info[ "status" ] == "N":
                i_img += 1
            else:
                img_ok = True

        new_img = np.zeros( [ conf.stim.img_diam_pix ] * 2 + [ 3 ] )

        # first set the coherent patches
        img_path = os.path.join( paths.img_db,
                                 coh_img_info[ "album" ],
                                 coh_img_info[ "image" ]
                               )

        coh_img = ns_patches.stimulus.load_img( img_path )

        for i_coh_patch in coh_patches:
            for i_chan in xrange( 3 ):
                new_img[ ..., i_chan ] += ( coh_img[ ..., i_chan ] *
                                            ( masks[ i_coh_patch, ... ] > 0 )
                                          )

        # now set the incoherent patches
        for i_incoh_patch in incoh_patches:

            valid_img_db = img_db_info[ img_db_info[ "status" ] == "Y" ]

            i_incoh_img = np.random.choice( len( valid_img_db ) )

            incoh_path = os.path.join( paths.img_db,
                                       valid_img_db[ i_incoh_img ][ "album" ],
                                       valid_img_db[ i_incoh_img ][ "image" ]
                                     )

            incoh_img = ns_patches.stimulus.load_img( incoh_path )

            for i_chan in xrange( 3 ):
                new_img[ ..., i_chan ] += ( incoh_img[ ..., i_chan ] *
                                            ( masks[ i_incoh_patch, ... ] > 0 )
                                          )

        stim.setTex( new_img )

        change_image = False
        mask_off = False

        while not change_image:

            map( lambda x : x.draw(), fix_stim + [ stim ] )

            stat_txt = [ "ID: {id:d}".format( id = coh_img_info[ "id" ] ),
                         "Album: " + coh_img_info[ "album" ],
                         "Image: " + coh_img_info[ "image" ],
                         "Status: " + coh_img_info[ "status" ],
                       ]

            status.setText( "\t".join( stat_txt ) )

            status.draw()

            win.flip()

            keys = psychopy.event.waitKeys()

            for key in keys:

                if key == "q":
                    img_db_info[ i_img ] = coh_img_info
                    keep_going = False
                    change_image = True

                elif key in [ "lshift", "rshift" ]:
                    img_db_info[ i_img ] = coh_img_info
                    i_img -= 1
                    change_image = True

                elif key in [ "space", "r" ]:
                    if not mask_off:
                        img_db_info[ i_img ] = coh_img_info
                        i_img += 1
                        change_image = True

                elif key == "y":
                    coh_img_info[ "status" ] = "Y"

                elif key == "n":
                    coh_img_info[ "status" ] = "N"

                elif key in [ "m", "b" ]:
                    if not mask_off:
                        temp_mask = stim._maskName.copy()
                        stim.setTex( coh_img )
                        stim.setMask( np.ones( coh_img.shape ) )
                        mask_off = True
                    else:
                        stim.setTex( new_img )
                        stim.setMask( temp_mask )
                        mask_off = False

                if change_image == True:
                    i_trial = np.mod( i_trial + 1, 80 )

    win.close()

    if save_changes:
        save_img_db_info( img_db_file = paths.img_db_info,
                          img_db_info = img_db_info
                        )


def read_img_db_info( img_db_file ):
    """Reads the image database file

    Parameters
    ----------
    img_db_file : string
        Path to the image database file

    Returns
    -------
    img_db_info : bmatch.prepare.prepare_images.write_img_db_info
        Image database info array

    """

    ( dtype, _ ) = get_img_db_dtype()

    img_db_info = np.loadtxt( fname = img_db_file,
                              dtype = dtype,
                              delimiter = "\t"
                            )

    return img_db_info


def get_img_db_dtype():
    """Returns the data type and formatting info for the image data.

    Returns
    -------
    dtype : numpy dtype
        The datatype info for each column in the data.
    file_formats : list of strings
        The formatting info to convert each column to text.

    """


    dtype = np.dtype( [ ( "id", "int" ),
                        ( "album", "|S255" ),
                        ( "image", "|S255" ),
                        ( "status", "|S1" ),
                      ]
                    )

    file_formats = [ "%i", "%s", "%s", "%s" ]

    return ( dtype, file_formats )


def save_img_db_info( img_db_file,
                      img_db_loc = None,
                      img_db_info = None,
                      album_list = None
                    ):
    """Saves the image database info to a file.

    Parameters
    ----------
    img_db_file : string
        Path to save the image database file.
    img_db_loc : string, optional
        Path to the image database (necessary if not passing `img_db_info`.
    img_db_info : img_db object or None, optional
        Image database information to save. If ``None``, an empty database is
        saved.
    album_list : list of strings, optional
        Only consider images from particular albums. Only useful in conjunction
        with `img_db_info` as `None`.

    """

    ( dtype, file_formats ) = get_img_db_dtype()

    if img_db_info is None:

        img_db_files = stimuli.mcgill_db.local_db_info( db_path = img_db_loc )

        ( albums, files ) = map( list, zip( *[ os.path.split( path ) for path in img_db_files ] ) )

        i_valid = np.where( [ album in album_list for album in albums ] )[ 0 ]

        img_db_files = [ img_db_files[ i_img ] for i_img in i_valid ]
        ( albums, files ) = map( list, zip( *[ os.path.split( path ) for path in img_db_files ] ) )

        img_db_info = np.empty( ( len( img_db_files ) ), dtype = dtype )

        img_db_info[ "id" ] = np.arange( len( img_db_files ) )
        img_db_info[ "album" ] = albums
        img_db_info[ "image" ] = files
        img_db_info[ "status" ] = "A"

    np.savetxt( fname = img_db_file,
                X = img_db_info,
                fmt = file_formats,
                header = "\t".join( dtype.names ),
                delimiter = "\t"
              )
