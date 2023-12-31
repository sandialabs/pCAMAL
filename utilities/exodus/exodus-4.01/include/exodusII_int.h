/*
//-------------------------------------------------------------------------
//
// Filename      : exodusII_int.h 
//
// Purpose       : File contains internal definitions for exodusII file format
//                 output needed within CUBIT to bypass certain very
//                 inefficient calls in the exodusII
//
// Special Notes : WARNING:: THIS IS A MODIFIED VERSION OF THE ORIGINAL
//                 EXODUS INCLUDE FILE, MODIFIED TO WORK IN C++ WITH THE
//                 CUBIT SOURCE CODE.
//
// Creator       : Ray J. Meyers
//
// Creation Date : 08/29/00
//
// Owner         : Ray J. Meyers
//-------------------------------------------------------------------------

// Using netcdf means we need to include the license file for it, so
// here it is:

//        Copyright 1993-2000 University
//         Corporation for Atmospheric
//         Research/Unidata

//         Portions of this software were developed by the Unidata
//         Program at the University Corporation for Atmospheric
//         Research.

//         Access and use of this software shall impose the following
//         obligations and understandings on the user. The user is
//         granted the right, without any fee or cost, to use, copy,
//         modify, alter, enhance and distribute this software, and
//         any derivative works thereof, and its supporting
//         documentation for any purpose whatsoever, provided that
//         this entire notice appears in all copies of the software,
//         derivative works and supporting documentation.  Further,
//         UCAR requests that the user credit UCAR/Unidata in any
//         publications that result from the use of this software or
//         in any product that includes this software. The names UCAR
//         and/or Unidata, however, may not be used in any advertising
//         or publicity to endorse or promote any products or
//         commercial entity unless specific written permission is
//         obtained from UCAR/Unidata. The user also understands that
//         UCAR/Unidata is not obligated to provide the user with any
//         support, consulting, training or assistance of any kind
//         with regard to the use, operation and performance of this
//         software nor to provide the user with any updates,
//         revisions, new versions or "bug fixes."

//         THIS SOFTWARE IS PROVIDED BY UCAR/UNIDATA "AS IS" AND ANY
//         EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//         THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//         PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
//         UCAR/UNIDATA BE LIABLE FOR ANY SPECIAL, INDIRECT OR
//         CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
//         FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
//         CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
//         OR IN CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS
//         SOFTWARE.
*/
/*****************************************************************************
*
* exodusII_int.h - ExodusII header file for internal Exodus call use only
*
* author - Sandia National Laboratories
*          Vic Yarberry    - Added headers and error logging
*
*          
* environment - UNIX
*
* revision history - 
*
*
****************************************************************************
*/

#ifndef EXODUS_II_INT_HDR
#define EXODUS_II_INT_HDR

#include "netcdf.h"
#ifdef NT
#include <stdlib.h>
#endif
#if defined __STDC__ || defined __cplusplus
#include <stdlib.h>
#endif

#include <stdio.h>

#if defined __STDC__ || defined __cplusplus
#define PROTO_ARGS(proto) proto
#define CONST_CHAR (const char *)
#define CONST const
#define VOID_PTR (void *)
#else
#define PROTO_ARGS(proto) ()
#define CONST_CHAR
#define VOID_PTR
#define CONST
#endif

#ifdef __cplusplus
extern "C" {
#endif  

/* these should be defined in ANSI C, and probably C++, but just in case ... */

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef NULL
#define NULL 0
#endif

/* EXODUS II version number */

/* ExodusII file version */
#define EX_VERS	3.01
/* ExodusII access library version */
#define EX_API_VERS 4.01

/* Exodus error return codes - function return values:                      */
#define EX_FATAL        -1      /* fatal error flag def                     */
#define EX_NOERR         0      /* no error flag def                        */
#define EX_WARN          1      /* warning flag def                         */

/*
 * This file contains defined constants that are used internally in the
 * EXODUS II API.
 *
 * The first group of constants refer to netCDF variables, attributes, or 
 * dimensions in which the EXODUS II data are stored.  Using the defined 
 * constants will allow the names of the netCDF entities to be changed easily 
 * in the future if needed.  The first three letters of the constant identify 
 * the netCDF entity as a variable (VAR), dimension (DIM), or attribute (ATT).
 *
 * NOTE: The entity name should not have any blanks in it.  Blanks are
 *       technically legal but some netcdf utilities (ncgen in particular)
 *       fail when they encounter a blank in a name.
 *
 *      DEFINED CONSTANT	ENTITY NAME	DATA STORED IN ENTITY
 */
#define ATT_FILE_TYPE		"type"		/* obsolete                  */
#define ATT_TITLE		"title"		/* the database title        */
#define ATT_API_VERSION		"api_version"	/* the EXODUS II api vers #   */
#define ATT_API_VERSION_BLANK	"api version"	/* the EXODUS II api vers #   */
						/*  used for db version 2.01 */
						/*  and earlier              */
#define ATT_VERSION		"version"	/* the EXODUS II file vers # */
#define ATT_FILESIZE		"file_size"	/* 1=large, 0=normal */
#define ATT_FLT_WORDSIZE	"floating_point_word_size"
						/* word size of floating     */
						/* point numbers in file     */
#define ATT_FLT_WORDSIZE_BLANK	"floating point word size"
						/* word size of floating     */
						/* point numbers in file     */
						/* used for db version 2.01  */
						/* and earlier               */
#define DIM_NUM_NODES		"num_nodes"	/* # of nodes                */
#define DIM_NUM_DIM		"num_dim"	/* # of dimensions; 2- or 3-d*/
#define DIM_NUM_ELEM		"num_elem"	/* # of elements             */
#define DIM_NUM_EL_BLK		"num_el_blk"	/* # of element blocks       */
#define VAR_COORD		"coord"		/* nodal coordinates         */
#define VAR_COORD_X             "coordx"        /* X-dimension coordinate    */
#define VAR_COORD_Y             "coordy"        /* Y-dimension coordinate    */
#define VAR_COORD_Z             "coordz"        /* Z-dimension coordinate    */
#define VAR_NAME_COOR		"coor_names"	/* names of coordinates      */
#define VAR_STAT_EL_BLK		"eb_status"	/* element block status      */
#define VAR_ID_EL_BLK		"eb_prop1"	/* element block ids props   */
#define ATT_NAME_ELB		"elem_type"	/* element type names for    */
						/*   each element block      */
#define DIM_NUM_EL_IN_BLK(num)	ex_catstr("num_el_in_blk",num)
						/* # of elements in element  */
 						/*   block num               */
#define DIM_NUM_NOD_PER_EL(num)	ex_catstr("num_nod_per_el",num)
						/* # of nodes per element in */
						/*   element block num       */
#define DIM_NUM_ATT_IN_BLK(num)	ex_catstr("num_att_in_blk",num)
						/* # of attributes in element*/
						/*   block num               */
#define VAR_CONN(num)		ex_catstr("connect",num)
						/* element connectivity for  */
						/*   element block num       */
#define VAR_ATTRIB(num)		ex_catstr("attrib",num)
						/* list of attributes for    */
						/*   element block num       */
#define VAR_EB_PROP(num)	ex_catstr("eb_prop",num)
						/* list of the numth property*/
						/*   for all element blocks  */
#define ATT_PROP_NAME		"name"		/* name attached to element  */
						/*   block, node set, side   */
						/*   set, element map, or    */
                                                /*   map properties          */
#define VAR_MAP			"elem_map"	/* element order map         */
                                                /* obsolete, replaced by     */
                                                /* VAR_ELEM_MAP(num)         */
#define DIM_NUM_SS		"num_side_sets"	/* # of side sets            */
#define VAR_SS_STAT		"ss_status"   	/* side set status           */

#ifdef CUBIT
#define VAR_SS_IDS		VAR_SS_PROP(1)	/* side set id properties    */
#else
#define VAR_SS_IDS		"ss_prop1"	/* side set id properties    */
#endif
#define DIM_NUM_SIDE_SS(num)	ex_catstr("num_side_ss",num)
                                                /* # of sides in side set num*/
#define DIM_NUM_DF_SS(num)	ex_catstr("num_df_ss",num)
                                               /* # of distribution factors */
                                               /* in side set num           */
/*#define DIM_NUM_NOD_SS(num)	ex_catstr("num_nod_ss",num) *** obsolete *** */
						/* # of nodes in side set num*/
#define VAR_FACT_SS(num)	ex_catstr("dist_fact_ss",num)
						/* the distribution factors  */
						/*   for each node in side   */
						/*   set num                 */
#define VAR_ELEM_SS(num)	ex_catstr("elem_ss",num)
						/* list of elements in side  */
						/*   set num                 */
#define VAR_SIDE_SS(num)        ex_catstr("side_ss",num)
                                                /* list of sides in side set */
#define VAR_SS_PROP(num)	ex_catstr("ss_prop",num)
						/* list of the numth property*/
						/*   for all side sets       */
#define DIM_NUM_NS		"num_node_sets"	/* # of node sets            */
#define DIM_NUM_NOD_NS(num)	ex_catstr("num_nod_ns",num)
						/* # of nodes in node set    */
						/*   num                     */
#define DIM_NUM_DF_NS(num)	ex_catstr("num_df_ns",num)
                                               /* # of distribution factors */
                                               /* in node set num           */
#define VAR_NS_STAT		"ns_status"	/* node set status           */

#ifdef CUBIT
#define VAR_NS_IDS		VAR_NS_PROP(1)	/* node set id properties    */
#else
#define VAR_NS_IDS		"ns_prop1"	/* node set id properties    */
#endif

#define VAR_NODE_NS(num)	ex_catstr("node_ns",num)
						/* list of nodes in node set */
						/*   num                     */
#define VAR_FACT_NS(num)	ex_catstr("dist_fact_ns",num)
						/* list of distribution      */
						/*   factors in node set num */
#define VAR_NS_PROP(num)	ex_catstr("ns_prop",num)
						/* list of the numth property*/
						/*   for all node sets       */
#define DIM_NUM_QA		"num_qa_rec"	/* # of QA records           */
#define VAR_QA_TITLE		"qa_records"	/* QA records                */
#define DIM_NUM_INFO		"num_info"	/* # of information records  */
#define VAR_INFO		"info_records"	/* information records       */
#define VAR_HIS_TIME		"time_hist"	/* obsolete                  */
#define VAR_WHOLE_TIME		"time_whole"	/* simulation times for whole*/
						/*   time steps              */
#define VAR_ELEM_TAB		"elem_var_tab"	/* element variable truth    */
						/*   table                   */
#define DIM_NUM_GLO_VAR		"num_glo_var"	/* # of global variables     */
#define VAR_NAME_GLO_VAR	"name_glo_var"	/* names of global variables */
#define VAR_GLO_VAR		"vals_glo_var"	/* values of global variables*/
#define DIM_NUM_NOD_VAR		"num_nod_var"	/* # of nodal variables      */
#define VAR_NAME_NOD_VAR	"name_nod_var"	/* names of nodal variables  */
#define VAR_NOD_VAR		"vals_nod_var"	/* values of nodal variables */
#define VAR_NOD_VAR_NEW(num)	ex_catstr("vals_nod_var",num)
                                                /* values of nodal variables */
#define DIM_NUM_ELE_VAR		"num_elem_var"	/* # of element variables    */
#define VAR_NAME_ELE_VAR	"name_elem_var"	/* names of element variables*/
#define VAR_ELEM_VAR(num1,num2)	ex_catstr2("vals_elem_var",num1,"eb",num2)
						/* values of element variable*/
						/*   num1 in element block   */
						/*   num2                    */
#define DIM_NUM_HIS_VAR		"num_his_var"	/* obsolete                  */
#define VAR_NAME_HIS_VAR	"name_his_var"	/* obsolete                  */
#define VAR_HIS_VAR		"vals_his_var"	/* obsolete                  */
#define DIM_STR			"len_string"	/* general dimension of      */
						/*   length MAX_STR_LENGTH   */
						/*   used for name lengths   */
#define DIM_LIN			"len_line"	/* general dimension of      */
 						/*   length MAX_LINE_LENGTH  */
						/*   used for long strings   */
#define DIM_N4			"four"		/* general dimension of      */
						/*   length 4                */
#define DIM_TIME		"time_step"	/* unlimited (expandable)    */
						/*   dimension for time steps*/
#define DIM_HTIME		"hist_time_step"/* obsolete                  */
#define VAR_ELEM_NUM_MAP        "elem_num_map"  /* element numbering map     */
                                                /* obsolete, replaced by     */
                                                /* VAR_ELEM_MAP(num)         */
#define VAR_NODE_NUM_MAP        "node_num_map"  /* node numbering map        */
                                                /* obsolete, replaced by     */
                                                /* VAR_NODE_MAP(num)         */
#define DIM_NUM_EM		"num_elem_maps"	/* # of element maps         */
#define VAR_ELEM_MAP(num)	ex_catstr("elem_map",num)
						/* the numth element map     */
#define VAR_EM_PROP(num)	ex_catstr("em_prop",num)
						/* list of the numth property*/
						/*   for all element maps    */
#define DIM_NUM_NM		"num_node_maps"	/* # of node maps            */
#define VAR_NODE_MAP(num)	ex_catstr("node_map",num)
						/* the numth node map        */
#define VAR_NM_PROP(num)	ex_catstr("nm_prop",num)
						/* list of the numth property*/
						/*   for all node maps       */

#define NUM_CFRAMES  "num_cframes"
#define NUM_CFRAME9  "num_cframes_9"
#define FRAME_COORDS "frame_coordinates"
#define FRAME_IDS    "frame_ids"
#define FRAME_TAGS   "frame_tags"


#define UNK			-1		/* unknown entity */
#ifndef CUBIT
#define TRIANGLE		1		/* Triangle entity */
#define QUAD			2		/* Quad entity */
#define HEX			3		/* Hex entity */
#define WEDGE			4		/* Wedge entity */
#define TETRA			5		/* Tetra entity */
#define TRUSS			6		/* Truss entity */
#define BEAM			7		/* Beam entity */
#define SHELL			8		/* Shell entity */
#define SPHERE			9		/* Sphere entity */
#define CIRCLE		       10		/* Circle entity */
#define TRISHELL               11               /* Triangular Shell */
#define PYRAMID                12               /* Pyramid */
#endif
/* Internal structure declarations */

struct list_item {              /* for use with ex_get_file_item */

  int exo_id;
  int value;
  struct list_item* next;
};

struct obj_stats {
  nclong *id_vals;
  nclong *stat_vals;
  long num;
  int exoid;
  int valid_ids;
  int valid_stat;
  struct obj_stats *next;
};


/* declare function prototypes.  recall that PROTO_ARGS() is a macro that
 * puts argument list in prototype if this is ANSI C.
 */

void ex_iqsort   PROTO_ARGS( (int v[], int iv[], int count ) );
char *ex_catstr PROTO_ARGS( (CONST char*, int) );
char *ex_catstr2 PROTO_ARGS( (CONST char*, int, CONST char*, int) );

enum convert_task { RTN_ADDRESS, 
                    READ_CONVERT, 
                    WRITE_CONVERT, 
                    WRITE_CONVERT_DOWN, 
                    WRITE_CONVERT_UP };

int ex_conv_ini PROTO_ARGS( ( int, int*, int*, int) );
void ex_conv_exit PROTO_ARGS( (int) );
nc_type nc_flt_code PROTO_ARGS( (int) );
int ex_comp_ws PROTO_ARGS( (int) );
void* ex_conv_array PROTO_ARGS( (int, int, const void*, int) );
int ex_get_cpu_ws();

void ex_rm_file_item_eb PROTO_ARGS( (int) );
void ex_rm_file_item_ns PROTO_ARGS( (int) );
void ex_rm_file_item_ss PROTO_ARGS( (int) );

extern struct list_item* eb_ctr_list;
extern struct list_item* ns_ctr_list;
extern struct list_item* ss_ctr_list;
extern struct list_item* em_ctr_list;
extern struct list_item* nm_ctr_list;

int ex_get_file_item PROTO_ARGS( ( int, struct list_item**) );
int ex_inc_file_item PROTO_ARGS( ( int,  struct list_item**) );
void ex_rm_file_item PROTO_ARGS( (int, struct list_item**) );

extern struct obj_stats* eb;
extern struct obj_stats* ns;
extern struct obj_stats* ss;
extern struct obj_stats* em;
extern struct obj_stats* nm;

int ex_get_side_set_node_list_len PROTO_ARGS( (int, int, int*) );
struct obj_stats *get_stat_ptr PROTO_ARGS( ( int, struct obj_stats**) );
void rm_stat_ptr PROTO_ARGS( (int, struct obj_stats**) );

int ex_id_lkup PROTO_ARGS( ( int exoid, char *id_type, int num) );
int ex_large_model PROTO_ARGS( (int exoid) );

#ifdef __cplusplus
}
#endif

#endif
