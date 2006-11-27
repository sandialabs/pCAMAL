/*****************************************************************************
*
* exodusII.h - Exodus II include file, for general use
*
* author - Sandia National Laboratories
*          Larry A. Schoof - Original
*          James A. Schutt - 8 byte float and standard C definitions
*          Vic Yarberry    - Added headers and error logging
*
*          
* environment - UNIX
*
* exit conditions - 
*
* revision history - 
*
*  Revision 1.1.1.1  2003/04/28 19:12:25  rakerr
*  importing exodusII files
*
*  Revision 1.24  2002/12/04 16:47:23  gdsjaar
*  Version changed to 3.26.
*  Garth Reese modifications:
*
*  Added support for coordinate frames: Coordinate frames are stored in
*  the database as a series of three points (defined in the basic
*  cartesian coordinate system). The first of these points describes the
*  origin of the new system. The second point lies on the 3 axis (or Z
*  axis) of the frame. The third point is in the 1-3 (xz) plane. Each
*  coordinate frame is identified by a unique, integer coordinate ID, and
*  by a character tag indicating whether the frame is rectangular
*  cartesian R , cylindrical C , or spherical S .
*
*  Two new api functions added:
*
*  Int ex_put_coordinate_frames(exoid, nframes, cf_ids, pt_coordinates, tags);
*
*  	int exoid   -- EXODUS file ID returned from a previous call to ex_create or ex_open.
*
*  	int nframes -- The number of coordinate frames to write
*
*  	const int cf_ids[] -- The (nframes) coordinate frame Ids. Integers greater than 0.
*
*  	const void* pt_coordinates The (9*nframes) coordinates of the
*  		three points defining each coordinate axis. The first three
*  		values are the origin of the first frame. The next three
*  		values are the coordinates of a point on the 3 rd axis of the
*  		first frame. The next three values are the coordinates of a
*  		point in the plane of the 1-3 axis. The pattern is repeated
*  		for each frame.
*
*  	const char* tags -- The (nframes) character tags associated with each coordinate frame.
*
*  Int ex_get_coordinate_frames(exoid, nframes, cf_ids, pt_coordinates, tags);
*
*  	int exoid --- EXODUS file ID returned from a previous call to ex_create or ex_open.
*
*  	int &nframes --- The number of coordinate frames to read
*
*  	const int cf_ids[] --- The (nframes) coordinate frame Ids. If
*  		cf_ids is NULL, no data will be returned in this or any other
*  		array. Only nframes will be modified. Otherwise, space must be
*  		allocated to store nframes integers before making this call.
*
*  	const void* pt_coordinates --- The (9*nframes) coordinates of the
*  		three points defining each coordinate axis. The first three
*  		values are the origin of the first frame. The next three
*  		values are the coordinates of a point on the 3 rd axis of the
*  		first frame. The next three values are the coordinates of a
*  		point in the plane of the 1-3 axis. The pattern is repeated
*  		for each frame. If cf_ids is null, no data will be returned in
*  		this array. Otherwise, space must be allocated for 9*nframes
*  		floating point values. The size of the allocation depends upon
*  		the compute word size.
*
*  	const char* tags --- The (nframes) character tags associated with
*  		each coordinate frame. If cf_ids is null, no data will be
*  		returned in this array. Otherwise, space must be allocated for
*  		nframes characters.
*
*  Revision 1.23  2001/09/20 21:51:53  gdsjaar
*  [Version changed to 3.19]
*  These modifications add a couple new routines that make it possible to
*  reduce the number of ncredef (netcdf redefinition) mode changes when
*  creating and writing a new exodusII file.  Switching into the ncredef
*  mode results in a complete rewrite of the file which can become very
*  expensive for large files.
*
*  Two new routines were added and a couple were modified.  The added
*  routines are:
*
*  int ex_put_concat_elem_block (int    exoid,
*  			      int*   elem_blk_id,
*  			      char *elem_type[],
*  			      int*   num_elem_this_blk,
*  			      int*   num_nodes_per_elem,
*  			      int*   num_attr,
*  			      int    define_maps)
*
*  	elem_blk_id 	   -- array of element block ids,
*  	elem_type	   -- array of element types,
*  	num_elem_this_blk  -- array of element block element counts,
*  	num_nodes_per_elem -- array of element node counts,
*  	num_attr	   -- array of element attribute counts,
*  	define_maps	   -- if != 0, define dimensions for element
*  			      and node numbering maps, else don't.
*
*  If this function is used, there will be a single ncredef call instead
*  of 1 per element block, 1 for each of the node/element numbering maps.
*
*  ex_put_concat_var_param
*
*  int ex_put_concat_var_param (int   exoid, int   num_g, int   num_n, int   num_e,
*  			     int   num_elem_blk, int  *elem_var_tab)
*
*  	int	exoid		exodus file id
*  	int     num_g   	global variable count
*  	int	num_n		nodal variable count
*  	int	num_e		element variable count
*  	int     num_elem_blk	number of element blocks
*  	int*    elem_var_tab	element variable truth table array
*
*  This function will result in 1 ncredef call instead of 1 per call to
*  ex_put_var_param and 1 in the call to ex_put_elem_var_tab.
*
*  The ex_put_concat_node_set function was modified.  If the
*  node_sets_node_index is passed as NULL (5th argument), then it only
*  defines the needed variables for the nodesets and does not write the
*  nodeset data; it can be written later with a call to ex_put_node_set.
*  This is useful if your code does not have all the data it would need
*  to define and write the nodesets at the same time, but want the
*  efficiency provided by the ex_put_concat_node_set routine.
*
*  A similar modification was made to ex_put_concat_side_set.  If the
*  side_sets_elem_index argument is passed as NULL (5th argument), then
*  it will only define the sidesets and not write the data.
*
*  This has been tested on several SEACAS and other codes and no problems
*  have been detected.
*
*  Revision 1.22  2000/09/06 15:55:50  gdsjaar
*  Modified exerrval and exoptval to be extern "C" for use with C++ codes.
*
*  Revision 1.21  1999/08/31 15:43:28  gdsjaar
*
*  Cleaned up some things in the exodusII cbind/include and cbind/src
*  direcotries to eliminate the warnings compiling with 'gcc -Wall' (full
*  warnings).
*
*  The changes are of the following form:
*
*  1. Missing or extra arguments to print/sprintf calls
*  2. Missing function prototypes in exodusII.h
*  +extern int ex_get_num_props             PROTO_ARGS( (int, int) );
*  (exodusII.h)
*
*  3. Missing function prototypes and include in exodusII_int.h:
*  +#include <stdio.h>
*  +void ex_iqsort   PROTO_ARGS( (int v[], int iv[], int left, int right)
*  );
*  +int ex_get_side_set_node_list_len PROTO_ARGS( (int, int, int*) );
*  +int ex_id_lkup PROTO_ARGS( ( int exoid, char *id_type, int num) );
*
*  4. Unused variables
*
*  Revision 1.20  1999/01/18 19:03:54  laschoo
*  changed EX_OK return code to EX_NOERR to prevent a conflict with defined
*  constant EX_OK in sysexits.h and unistd.h on some systems (i.e., SGIs)
*
*  Revision 1.19  1998/03/31 17:21:48  laschoo
*  changes for netCDF 3.4
*    added MAX_HEADER_SIZE constant;  the ex_put_init () function will not
*    allocate more than this for the header
*
*    modified EXODUS error codes to conform to (new) netCDF 3.4 error conventions;
*    positive errors are bad (fatal); negative are informative (warnings);
*    application codes are isolated from this behavior since the EXODUS functions
*    still return EX_OK (=0) for success, EX_FATAL (=-1) for failure, and
*    EX_WARN (=1) for success with warning;  the EXODUS error codes are printed
*    by the ex_err () function
*
*  Revision 1.18  1997/05/13 14:02:55  laschoo
*  added function prototype for ex_get_side_set_node_list
*
*  Revision 1.17  1996/12/24 19:46:03  laschoo
*  modified to allow multiple node and element maps;
*  changed API version to 2.09 and file version to 2.03
*
*  Revision 1.16  1996/08/12 16:24:05  laschoo
*  modified itol and ltoi function prototypes to use nclong for netcdf 2.4.2
*
*  Revision 1.15  1996/07/09 22:03:03  laschoo
*  changed version to 2.05
*
*  Revision 1.14  1995/09/20 17:37:31  mksmith
*  Upgrade to version 2.03
*
 * Revision 1.13  1994/08/24  15:29:49  laschoo
 * new API version 2.02 which includes Alpha support
 *
 * Revision 1.12  1994/03/30  15:37:42  vryarbe
 * Updated for V2.01
 *
 * Revision 1.11  1993/11/18  18:54:16  vryarbe
 * changed name of options flag to exoptval
 *
 * Revision 1.10  1993/10/18  19:46:34  vryarbe
 * removed const declaration from passed strings
 *
 * Revision 1.9  1993/09/24  21:08:45  vryarbe
 * added definitions for new inquire parameters
 *
 * Revision 1.8  1993/09/13  20:39:45  vryarbe
 * added new parameters for inquiry and a new error code
 *
 * Revision 1.6  1993/08/30  16:23:29  vryarbe
 * added return codes EX_WARN and EX_OK
 *
 * Revision 1.5  1993/08/27  22:52:52  vryarbe
 * modfied for use with property functions
 *
 * Revision 1.4  1993/08/23  20:44:33  vryarbe
 * Added CONST and some new error msg definitions.
 *
 * Revision 1.3  1993/07/08  21:49:06  vryarbe
 * bug fixes to date
 *
 * Revision 1.2  1993/07/01  22:27:41  vryarbe
 * updated header
 *
*
*****************************************************************************/
#include "netcdf.h"
#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0 
#endif

#ifndef EXODUS_II_HDR
#define EXODUS_II_HDR

/* NOTE:
 *	as of 21 August 1992, the EXODUS II library, C binding version, is
 *	compiled under C only.  However, this header file is designed so that
 *	it can be included into, and the library linked with, a C++ program.
 */

#if defined __STDC__ || defined __cplusplus
#define PROTO_ARGS(proto) proto
#define CONST_CHAR (const char *)
#define VOID_PTR (void *)
#else
#define PROTO_ARGS(proto) ()
#define CONST_CHAR
#define VOID_PTR
#endif

/* need following extern if this include file is used in a C++ program, to
 * keep the C++ compiler from mangling the function names.
 */
#ifdef __cplusplus
extern "C" {
#endif

/*
 * The following are miscellaneous constants used in the EXODUS II API.
 */

#define EX_NOCLOBBER	0
#define EX_CLOBBER	1
#define EX_NORMAL_MODEL 	2 /* disable mods that permit storage of larger models */
#define EX_LARGE_MODEL  	4 /* enable mods that permit storage of larger models */

#define EX_READ		0
#define EX_WRITE	1

#define EX_INQ_FILE_TYPE	1		/* inquire EXODUS II file type*/
#define EX_INQ_API_VERS		2		/* inquire API version number */
#define EX_INQ_DB_VERS		3		/* inquire database version   */
						/*   number                   */
#define EX_INQ_TITLE		4		/* inquire database title     */
#define EX_INQ_DIM		5		/* inquire number of          */
						/*   dimensions               */
#define EX_INQ_NODES		6		/* inquire number of nodes    */
#define EX_INQ_ELEM		7               /* inquire number of elements */
#define EX_INQ_ELEM_BLK		8		/* inquire number of element  */
						/*   blocks                   */
#define EX_INQ_NODE_SETS	9		/* inquire number of node sets*/
#define EX_INQ_NS_NODE_LEN	10		/* inquire length of node set */
						/*   node list                */
#define EX_INQ_SIDE_SETS	11		/* inquire number of side sets*/
#define EX_INQ_SS_NODE_LEN	12		/* inquire length of side set */
						/*   node list                */
#define EX_INQ_SS_ELEM_LEN	13		/* inquire length of side set */
						/*   element list             */
#define EX_INQ_QA		14		/* inquire number of QA       */
						/*   records                  */
#define EX_INQ_INFO		15		/* inquire number of info     */
						/*   records                  */
#define EX_INQ_TIME		16		/* inquire number of time     */
						/*   steps in the database    */
#define EX_INQ_EB_PROP          17              /* inquire number of element  */
                                                /*   block properties         */
#define EX_INQ_NS_PROP          18              /* inquire number of node set */
                                                /*   properties               */
#define EX_INQ_SS_PROP          19              /* inquire number of side set */
#define EX_INQ_NS_DF_LEN	20		/* inquire length of node set */
						/*   distribution factor  list*/
#define EX_INQ_SS_DF_LEN	21		/* inquire length of node set */
						/*   distribution factor  list*/
#define EX_INQ_LIB_VERS		22		/* inquire API Lib vers number*/
#define EX_INQ_EM_PROP          23              /* inquire number of element  */
                                                /*   map properties           */
#define EX_INQ_NM_PROP          24              /* inquire number of node     */
                                                /*   map properties           */
#define EX_INQ_ELEM_MAP		25		/* inquire number of element  */
						/*   maps                     */
#define EX_INQ_NODE_MAP		26		/* inquire number of node     */
						/*   maps                     */

/*   properties               */
#define EX_ELEM_BLOCK		1		/* element block property code*/
#define EX_NODE_SET		2		/* node set property code     */
#define EX_SIDE_SET		3		/* side set property code     */
#define EX_ELEM_MAP		4		/* element map property code  */
#define EX_NODE_MAP		5		/* node map property code     */

/*   max string lengths; constants that are used as netcdf dimensions must be
     of type long       */
#define MAX_STR_LENGTH		32L
#define MAX_VAR_NAME_LENGTH	20
#define MAX_LINE_LENGTH		80L
#define MAX_ERR_LENGTH          256

/*   for netCDF 3.4, we estimate the size of the header; 
     if estimate is larger than this max, set the estimate to this max;
     I've never measured a header larger than 20K   */
#define MAX_HEADER_SIZE         30000

/* declare function prototypes.  recall that PROTO_ARGS() is a macro that
 * puts argument list in prototype if this is ANSI C or C++.
 */

/* routines for file initialization i/o */

extern int ex_create			PROTO_ARGS( (const char*,
						     int, int*, int*) );
extern int ex_open			PROTO_ARGS( (const char*,
						     int, int*, int*, float*) );
extern int ex_close			PROTO_ARGS( (int) );
extern void ex_err			PROTO_ARGS( (const char*, const char*, int) );
extern void ex_opts			PROTO_ARGS( (int) );
extern int ex_update			PROTO_ARGS( (int) );

extern int ex_put_init			PROTO_ARGS( (int, const char*, int, int,
						     int, int, int, int) );
extern int ex_get_init			PROTO_ARGS( (int, char*, int*, int*,
						     int*, int*, int*, int*) );

extern int ex_put_qa			PROTO_ARGS( (int,int, char*[][4]));
extern int ex_get_qa			PROTO_ARGS( (int, char*[][4]) );

extern int ex_put_info			PROTO_ARGS( (int, int, char*[]) );
extern int ex_get_info			PROTO_ARGS( (int, char*[]) );

/* routines for model description i/o */

extern int ex_put_coord			PROTO_ARGS( (int, const void*, const void*,
						     const void*) );
extern int ex_get_coord			PROTO_ARGS( (int, void*, void*,
						     void*) );

extern int ex_put_coord_names		PROTO_ARGS( (int, char*[]) );
extern int ex_get_coord_names		PROTO_ARGS( (int, char*[]) );

extern int ex_put_map			PROTO_ARGS( (int, const int*) );
extern int ex_get_map			PROTO_ARGS( (int, int*) );

extern int ex_put_elem_block		PROTO_ARGS( (int, int, const char*, int,
						     int, int) );
extern int ex_get_elem_block		PROTO_ARGS( (int, int, char*, int*,
						     int*, int*) );
extern int ex_put_concat_elem_block	PROTO_ARGS( (int, int*, char*[], int*,
						     int*, int*, int) );

extern int ex_get_elem_blk_ids		PROTO_ARGS( (int, int*) );

extern int ex_put_elem_conn		PROTO_ARGS( (int, int, const int*) );
extern int ex_get_elem_conn		PROTO_ARGS( (int, int, int*) );

extern int ex_put_elem_attr		PROTO_ARGS( (int, int, const void*) );
extern int ex_get_elem_attr		PROTO_ARGS( (int, int, void*) );

extern int ex_put_node_set_param	PROTO_ARGS( (int, int, int, int) );
extern int ex_get_node_set_param	PROTO_ARGS( (int, int, int*, int*) );

extern int ex_put_node_set		PROTO_ARGS( (int, int, const int*) );
extern int ex_get_node_set		PROTO_ARGS( (int, int, int*) );

extern int ex_put_node_set_dist_fact	PROTO_ARGS( (int, int, const void*) );
extern int ex_get_node_set_dist_fact	PROTO_ARGS( (int, int, void*) );

extern int ex_get_node_set_ids		PROTO_ARGS( (int, int*) );

extern int ex_put_concat_node_sets	PROTO_ARGS( (int, int*, int*, int*,
						     int*, int*, int*, void*) );
extern int ex_get_concat_node_sets      PROTO_ARGS( (int, int*, int*, int*,
						     int*, int*, int*, void*) );

extern int ex_put_side_set_param	PROTO_ARGS( (int, int, int, int) );
extern int ex_get_side_set_param	PROTO_ARGS( (int, int, int*, int*) );

extern int ex_put_side_set		PROTO_ARGS( (int, int, const int*, const int*) );
extern int ex_get_side_set		PROTO_ARGS( (int, int, int*, int*) );
extern int ex_put_side_set_dist_fact	PROTO_ARGS( (int, int, const void*) );
extern int ex_get_side_set_dist_fact	PROTO_ARGS( (int, int, void*) );
extern int ex_get_side_set_ids		PROTO_ARGS( (int, int*) );
extern int ex_get_side_set_node_list	PROTO_ARGS( (int, int, int*, int*) );

extern int ex_put_prop_names		PROTO_ARGS( (int, int, int, char**) );
extern int ex_get_prop_names		PROTO_ARGS( (int, int, char**) );

extern int ex_put_prop			PROTO_ARGS( (int, int, int, const char*,
						     int) );
extern int ex_get_prop			PROTO_ARGS( (int, int, int, const char*,
						     int*) );

extern int ex_put_prop_array		PROTO_ARGS( (int, int, const char*, const int*) );
extern int ex_get_prop_array		PROTO_ARGS( (int, int, const char*, int*) );

extern int ex_put_concat_side_sets	PROTO_ARGS( (int, const int*, const int*, const int*,
						     const int*, const int*, const int*, const int*,
                                                     const void* ) );
extern int ex_get_concat_side_sets	PROTO_ARGS( (int, int*, int*, int*,
						     int*, int*, int*, int*,
                                                     void* ) );
extern int ex_cvt_nodes_to_sides	PROTO_ARGS( (int, int*, int*, int*,
						     int*, int*, int*, int*) );
extern int ex_put_coordinate_frames     PROTO_ARGS( (int, int, const int*, 
						    const void*, const char*));

extern int ex_get_coordinate_frames     PROTO_ARGS( (int, int*, int*, void*,
						     char*) );

/* routines for analysis results i/o */

extern int ex_put_var_param		PROTO_ARGS( (int, const char*, int) );
extern int ex_get_var_param		PROTO_ARGS( (int, const char*, int*) );

extern int ex_put_concat_var_param      PROTO_ARGS( (int, int, int, int, int, int*) );
						      
extern int ex_put_var_names		PROTO_ARGS( (int, const char*, int,
						     char*[]) );
extern int ex_get_var_names		PROTO_ARGS( (int, const char*, int,
						     char*[]) );

extern int ex_put_var_name		PROTO_ARGS( (int, const char*, int,
						     const char*) );
extern int ex_get_var_name		PROTO_ARGS( (int, const char*, int,
						     char*) );

extern int ex_put_elem_var_tab		PROTO_ARGS( (int, int, int, int*) );
extern int ex_get_elem_var_tab		PROTO_ARGS( (int, int, int, int*) );

extern int ex_put_glob_vars		PROTO_ARGS( (int, int, int, const void*) );
extern int ex_get_glob_vars		PROTO_ARGS( (int, int, int, void*) );

extern int ex_get_glob_var_time		PROTO_ARGS( (int, int, int, int,
						     void*) );

extern int ex_put_nodal_var		PROTO_ARGS( (int, int, int, int,
						     const void*) );
extern int ex_get_nodal_var		PROTO_ARGS( (int, int, int, int,
						     void*) );

extern int ex_get_nodal_var_time	PROTO_ARGS( (int, int, int, int, int,
						     void*) );

extern int ex_put_elem_var		PROTO_ARGS( (int, int, int, int, int,
						     const void*) );
extern int ex_get_elem_var		PROTO_ARGS( (int, int, int, int, int,
						     void*) );

extern int ex_get_elem_var_time		PROTO_ARGS( (int, int, int, int, int,
						     void*) );

extern int ex_put_time			PROTO_ARGS( (int, int, const void*) );
extern int ex_get_time			PROTO_ARGS( (int, int, void*) );

extern int ex_get_all_times		PROTO_ARGS( (int, void*) );

extern int ex_inquire			PROTO_ARGS( (int, int, int*, void*,
						     char*) );
extern int ex_get_num_props             PROTO_ARGS( (int, int) );

extern int ex_put_elem_num_map          PROTO_ARGS( (int, const int*) );
extern int ex_get_elem_num_map          PROTO_ARGS( (int, int*) );

extern int ex_put_node_num_map          PROTO_ARGS( (int, const int*) );
extern int ex_get_node_num_map          PROTO_ARGS( (int, int*) );

extern int ex_put_map_param             PROTO_ARGS( (int, int, int) );
extern int ex_get_map_param             PROTO_ARGS( (int, int*, int*) );

extern int ex_put_elem_map              PROTO_ARGS( (int, int, const int*) );
extern int ex_get_elem_map              PROTO_ARGS( (int, int, int*) );

extern int ex_put_node_map              PROTO_ARGS( (int, int, const int*) );
extern int ex_get_node_map              PROTO_ARGS( (int, int, int*) );

extern nclong *itol                     PROTO_ARGS( (const int*, int) ); 

extern int ltoi                         PROTO_ARGS( (const nclong*, int*, int) );

extern int ex_copy                     PROTO_ARGS( (int, int) );

extern int cpy_att                     PROTO_ARGS( (int, int, int, int) );

extern int cpy_var_def                  PROTO_ARGS( (int, int, int, char*) );

extern int cpy_var_val                 PROTO_ARGS( (int, int, char*) );

/* ERROR CODE DEFINITIONS AND STORAGE                                       */
extern int exerrval;            /* shared error return value                */
extern int exoptval;            /* error reporting flag (default is quiet)  */

#ifdef __cplusplus
}				/* close brackets on extern "C" declaration */
#endif

#endif

/* ex_opts function codes - codes are OR'ed into exopts                     */
#define EX_VERBOSE      1	/* verbose mode message flag                */
#define EX_DEBUG        2       /* debug mode def                           */
#define EX_ABORT        4       /* abort mode flag def                      */

/* Exodus error return codes - exerrval return values:                      */
#define EX_MEMFAIL	1000    /* memory allocation failure flag def       */
#define EX_BADFILEMODE	1001    /* bad file mode def                        */
#define EX_BADFILEID	1002    /* bad file id def                          */
#define EX_WRONGFILETYPE 1003   /* wrong file type for function             */
#define EX_LOOKUPFAIL	1004    /* id table lookup failed                   */
#define EX_BADPARAM	1005    /* bad parameter passed                     */
#define EX_NULLENTITY	-1006   /* null entity found                        */
#define EX_MSG		-1000   /* message print code - no error implied    */
#define EX_PRTLASTMSG   -1001   /* print last error message msg code        */

