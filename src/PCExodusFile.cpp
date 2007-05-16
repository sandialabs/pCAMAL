// file: PCExodusFile.cpp
// description: implementation of pCAMAL Exodus file interface
// author: Michael Stephenson

#include <set>
#include <string.h>
#include "PCExodusFile.hpp"
#include "PCSweepVolume.hpp"
#include "PCMLSweeper.hpp"

PCExodusFile::PCExodusFile(char* const file_name, pce::FileOp op)
        : exoID(0), cpuWord(8), fileWord(0), numDim(0), numNodes(0),
          numElems(0), numElemBlks(0), numNodeSets(0), numSideSets(0),
          mVersion(0.0)
{
  switch(op) {
      // open exising file for read only
    case pce::read:
        exoID = ex_open(file_name, EX_READ, &cpuWord, &fileWord,
                        &mVersion);
        read_init();
        break;

          // create a new file
    case pce::create:
        fileWord = 8;
        exoID = ex_create(file_name, EX_CLOBBER, &cpuWord, &fileWord);
        break;

          // open existing file for read/write
    case pce::write:
        exoID = ex_open(file_name, EX_WRITE, &cpuWord, &fileWord,
                        &mVersion);
        read_init();
        break;
    default:
        break;
  }

  if (exoID > 0) {
    strncpy(fileName, file_name, MAX_STR_LENGTH);
    fileName[MAX_STR_LENGTH] = '\0';
  }
  else {
    fileName[0] = '\0';
  }
}

PCExodusFile::~PCExodusFile()
{
  delete_sweep_volumes();
}
  
int PCExodusFile::get_num_hexes(int num_blks, int* num_hexes)
{
  if (exoID == 0 || num_blks < sweepVols.size())
    return 0;

  int total = 0;
  int i;
  for (i = 0; i < sweepVols.size(); i++) {
    num_hexes[i] = sweepVols[i]->get_num_hexes();
    total += num_hexes[i];
  }
  return total;
}

int PCExodusFile::get_num_sweep_vols()
{
  if (exoID == 0)
    return 0;

  return sweepVols.size();
}

void PCExodusFile::get_param(int& num_blks, int& num_node_sets, 
                             int& num_side_sets)
{
  num_blks = num_node_sets = num_side_sets = 0;
  if (exoID == 0)
    return;

  num_blks      = numElemBlks;
  num_node_sets = numNodeSets;
  num_side_sets = numSideSets;
}

void PCExodusFile::read_sweep_prop(int vol_id, int &block_id, int &sweep_id, 
                                   int &num_quads, int &nodes_per_hex)
{
  if (exoID == 0 || sweepVols.empty() ||
      vol_id < 0 || vol_id >= sweepVols.size()) {
    sweep_id  = 0;
    num_quads = 0;
    block_id  = 0;
    nodes_per_hex = 8;
  }
  else {
    sweep_id  = sweepVols[vol_id]->get_sweep_id();
    num_quads = sweepVols[vol_id]->get_num_quads();
    block_id  = sweepVols[vol_id]->get_elem_block_id();
    nodes_per_hex = sweepVols[vol_id]->get_nodes_per_hex();
  }  
}

bool PCExodusFile::get_node_sets(int num_points, int* node_ids, 
                                 int& new_node_sets, int* ns_id_array, 
                                 int* ns_cnts_array, int* ns_df_cnts_array, 
                                 int* ns_ptrs, int* ns_df_ptrs,
                                 int* &ns_list, double* &ns_df_list)
{
  bool debug_flag = false;
  
  if (exoID == 0 || numNodeSets == 0)
    return false;

    // allocate memory for lists
  int ns_len = 0;
  int ns_df_len = 0;
  float fdum;
  char cdum;
  int error = ex_inquire(exoID, EX_INQ_NS_NODE_LEN, &ns_len, &fdum, &cdum);
  if (error == 0) {
    ns_list = new int[ns_len];
    if (ns_list == NULL)
      return false;

    error = ex_inquire(exoID, EX_INQ_NS_DF_LEN, &ns_df_len, &fdum, &cdum);
  }

  if (error == 0 && ns_df_len > 0) {
    ns_df_list = new double[ns_df_len];
    if (ns_df_list == NULL) {
      delete [] ns_list;
      ns_list = NULL;
      return false;
    }
    
    error = ex_get_concat_node_sets(exoID, ns_id_array, ns_cnts_array,
                                    ns_df_cnts_array, ns_ptrs, ns_df_ptrs,
                                    ns_list, ns_df_list);
  }

  if (debug_flag && error == 0) {
    printf("Before conversion\n");
    print_concat_node_sets(ns_id_array, ns_cnts_array, ns_df_cnts_array,
                           ns_ptrs, ns_df_ptrs, ns_list, ns_df_list);

    double* xx = new double[numNodes];
    double* yy = new double[numNodes];
    double* zz = new double[numNodes];
    if (xx != NULL && yy != NULL && zz != NULL) {
      if (ex_get_coord(exoID, xx, yy, zz) == 0) {
        print_node_sets_coords(ns_id_array, ns_cnts_array, ns_ptrs, ns_list,
                               xx, yy, zz);
      }
    }
    
    delete [] zz;
    delete [] yy;
    delete [] xx;
  }

  if (error == 0) {
    int i, j;
    std::map<int,int> node_map;
    for (i = 0; i < num_points; i++) {
      node_map[node_ids[i]] = i + 1;
    }

      // convert to new node numbers
    int new_ns_len = 0;
    for (i = 0; i < numNodeSets; i++) {
      int begin = ns_ptrs[i];
      int end = begin + ns_cnts_array[i];
      for (j = begin; j < end; j++) {
        if (node_map.find(ns_list[j]) != node_map.end()) {
          ns_list[new_ns_len] = node_map[ns_list[j]];
          if (ns_df_len > 0)
            ns_df_list[new_ns_len] = ns_df_list[j];
          ++new_ns_len;
        }
      }
      ns_ptrs[i] = new_ns_len;
    }

    ns_cnts_array[0] = ns_ptrs[0];
    for (i = 1; i < numNodeSets; i++) {
      ns_cnts_array[i] = ns_ptrs[i] - ns_ptrs[i-1];
    }
    ns_ptrs[0] = 0;
    for (i = 1; i < numNodeSets; i++) {
      ns_ptrs[i] = ns_ptrs[i-1] + ns_cnts_array[i];
    }
    if (ns_df_len > 0) {
      for (i = 0; i < numNodeSets; i++) {
        ns_df_cnts_array[i] = ns_cnts_array[i];
        ns_df_ptrs[i] = ns_ptrs[i];
      }
    }
  
    if (debug_flag) {
      printf("After conversion\n");
      print_concat_node_sets(ns_id_array, ns_cnts_array, ns_df_cnts_array,
                             ns_ptrs, ns_df_ptrs, ns_list, ns_df_list);
    }

    if (new_ns_len == ns_len) {
      printf("No change to side set\n");
    }
  }

  return error == 0 ? true : false;
}

void PCExodusFile::print_concat_node_sets(int* ns_id_array, int* ns_cnts_array,
                                          int* ns_df_cnts_array, int* ns_ptrs,
                                          int* ns_df_ptrs, int* ns_list, 
                                          double* ns_df_list)
{
    // header info
  int i, j, k;
  for (i = 0; i < numNodeSets; i++) {
    printf("node set %d\n", i+1);
    printf("\tnode_set id     = %d\n", ns_id_array[i]);
    printf("\tnum nodes       = %d\n", ns_cnts_array[i]);
    printf("\tnum dist_fact   = %d\n", ns_df_cnts_array[i]);
    printf("\tindex node_set  = %d\n", ns_ptrs[i]);
    printf("\tindex dist_fact = %d\n", ns_df_ptrs[i]);
  }

  int ns_len = 0;
  int ns_df_len = 0;
  for (i = 0; i < numNodeSets; i++) {
    ns_len += ns_cnts_array[i];
    ns_df_len += ns_df_cnts_array[i];
  }

    // node set
  if (ns_len > 0) {
    for (i = 0; i <numNodeSets; i++) {
      printf("\nnode_ns%d\n", i+1);
      int begin = ns_ptrs[i];
      int end   = begin + ns_cnts_array[i];
      for (j = begin; j < end; j += 10) {
        printf("%4d:", j+1);
        for (k = 0; k < 10; k++) {
          if (j + k < end)
            printf(" %4d", ns_list[j+k]);
          else
            break;
        }
        printf("\n");
      }
    }
  }

    // distribution factors
  if (ns_df_len > 0) {
    for (i = 0; i <numNodeSets; i++) {
      printf("\ndist_fact_ns%d\n", i+1);
      int begin = ns_df_ptrs[i];
      int end   = begin + ns_df_cnts_array[i];
      for (j = begin; j < end; j += 5) {
        printf("%d:", j+1);
        for (k = 0; k < 5; k++) {
          if (j + k < end)
            printf(" %f", ns_df_list[j+k]);
          else
            break;
        }
        printf("\n");
      }
    }
  }
}

void PCExodusFile::print_node_sets_coords(int* ns_id_array, int* ns_cnts_array,
                                          int* ns_ptrs, int* ns_list,
                                          double* xx, double* yy, double* zz)
{
  int i, j;
  for (i = 0; i < numNodeSets; i++) {
    printf("%d: node_set%d\n", i+1, ns_id_array[i]);
    int begin = ns_ptrs[i];
    int end = begin + ns_cnts_array[i];
    for (j = begin; j < end; j++) {
      int k = ns_list[j] - 1;
      printf("%4d - %4d: %f %f %f\n", j+1, ns_list[j], xx[k], yy[k], zz[k]);
    }
  }
}


void PCExodusFile::read_sweep_coord(int vol_id, int& num_points, 
                                    double* &x_coor, double* &y_coor, 
                                    double* &z_coor, int* &node_ids)
{
  if (vol_id < 0 || vol_id > sweepVols.size()) {
    num_points = 0;
    return;
  }
  
    // retrieve all nodes
  double* xx = new double[numNodes];
  double* yy = new double[numNodes];
  double* zz = new double[numNodes];
  if (xx == NULL || yy == NULL || zz == NULL) {
    delete [] zz;
    delete [] yy;
    delete [] xx;
    return;
  }
  
  int error = ex_get_coord(exoID, xx, yy, zz);      

    // copy only coordinates in element blocks
  int* mark = new int[numNodes];
  if (mark == NULL) {
    delete [] mark;
    delete [] zz;
    delete [] yy;
    delete [] xx;
    return;
  }
  if (error == 0) {
      // get element blocks
    PCSweepVolume* vol = sweepVols[vol_id];
    std::vector<int> blk_ids;
    vol->get_ordered_ids(blk_ids);

      // mark all nodes in element blocks
    memset(mark, 0, numNodes * sizeof(int));
    
    int i;
    for (i = 0; i < blk_ids.size() && error == 0; i++) {
      int num_quads = vol->get_num_surf_quads(i);
      int conn[num_quads * 4];
      
      error = ex_get_elem_conn(exoID, blk_ids[i], conn);
      if (error == 0) {
        int j;
        for (j = 0; j < num_quads * 4; j++)
          mark[conn[j] - 1] = 1;
      }
    }
  }
    
    // count number of marked nodes
  num_points = 0;
  if (error == 0) {
    int i;
    for (i = 0; i < numNodes; i++) {
      if (mark[i] > 0)
        ++num_points;
    }

      // copy marked coordinates
    x_coor = new double[num_points];
    y_coor = new double[num_points];
    z_coor = new double[num_points];
    node_ids = new int[num_points];
    if (x_coor == NULL || y_coor == NULL || z_coor == NULL || 
        node_ids == NULL) {
      delete [] node_ids;
      delete [] z_coor;
      delete [] y_coor;
      delete [] x_coor;
    }

    else {
      int j = 0;
      for (i = 0; i < numNodes; i++) {
        if (mark[i] > 0) {
          x_coor[j] = xx[i];
          y_coor[j] = yy[i];
          z_coor[j] = zz[i];
          node_ids[j++] = i + 1;
//         printf("%d:%d %f %f %f\n", j, node_ids[j], xx[i], yy[i], zz[i]);
        }
      }
    }
  }

  delete [] mark;
  delete [] zz;
  delete [] yy;
  delete [] xx;
}

void PCExodusFile::read_sweep_conn(int vol_id, int num_points, int num_quads,
                                   int *node_ids, int *connect)
{
    // no Exodus file or no sweep volumes or index out of range
  if (exoID == 0 || sweepVols.empty() ||
      vol_id < 0 || vol_id >= sweepVols.size())
    return;

    // read the connectivity
  if (num_quads > 0 && node_ids != NULL && connect != NULL) {
    std::vector<int> blk_ids;
    int num_blks = sweepVols[vol_id]->get_ordered_ids(blk_ids);
    
    int* tc = connect;
    int error = 0;
    int j;
    for (j = 0; j < num_blks  && error == 0; j++) {
      int blk_id = blk_ids[j];
      error = ex_get_elem_conn(exoID, blk_id, tc);
      int num_blk_quads = sweepVols[vol_id]->get_num_surf_quads(j);
      tc += (num_blk_quads * 4);
    }
    
      // map node numbers for node_ids to sequential ids
    if (error == 0) {
      std::map<int,int> node_map;
      int j;
      for (j = 0; j < num_points; j++) {
        node_map[node_ids[j]] = j;
      }
      
      for (j = 0; j < num_quads * 4; j++) {
        int id = connect[j];
        connect[j] = node_map[id];
//         printf("%5d: %5d -> %5d\n", j, id, connect[j]);
      }
    }
  }
}

int PCExodusFile::read_sweep_surf_prop(int vol_id, int &num_src_surf,
                                        int &num_lnk_surf,
                                        int &num_tgt_surf)
{
  if (exoID == 0 || sweepVols.empty() ||
      vol_id < 0 || vol_id >= sweepVols.size()) {
    num_src_surf = 0;
    num_lnk_surf = 0;
    num_tgt_surf = 0;
  }
  else {
    int total = sweepVols[vol_id]->get_num_surfs(num_src_surf, num_lnk_surf,
                                                 num_tgt_surf);
  }
}

void PCExodusFile::read_sweep_surf_size(int vol_id, int num_surfs,
                                        int *num_surf_quads)
{
  if (exoID == 0 || sweepVols.empty() ||
      vol_id < 0 || vol_id >= sweepVols.size())
    return;
  
  if (num_surfs > 0 && num_surf_quads != NULL) {
    sweepVols[vol_id]->get_num_surf_quads(num_surfs, num_surf_quads);
  }
}

int PCExodusFile::put_param(int num_points_out, int num_hexes,
                            int num_node_sets, int num_side_sets)
{
  if (exoID == 0 || num_points_out == 0 || num_hexes == 0)
    return 1;

    // date and time info
  time_t result = time(NULL);
  struct tm *timeptr = localtime(&result);
  char camal_date[MAX_STR_LENGTH+1];
  strftime(camal_date, MAX_STR_LENGTH, "%D", timeptr);
  char camal_time[MAX_STR_LENGTH+1];
  strftime(camal_time, MAX_STR_LENGTH, "%H:%M:%S", timeptr);
  
    // write file parameters
  char title[MAX_LINE_LENGTH+1];
  snprintf(title, MAX_LINE_LENGTH, "pCAMAL(%s): %s: %s", 
           fileName, camal_date, camal_time);
  numDim      = 3;
  numNodes    = num_points_out;
  numElems    = num_hexes;
  numElemBlks = 1;
  numNodeSets = num_node_sets;
  numSideSets = num_side_sets;
  int error = ex_put_init(exoID, title, numDim, numNodes, numElems,
                          numElemBlks, numNodeSets, numSideSets);

    // write qa record
  if (error == 0) {
    char *qa[1][4];
    int i;
    for (i = 0; i < 4; i++)
      qa[0][i] = new char[MAX_STR_LENGTH+1];
    
      // write pCamal QA record
    strcpy(qa[0][0], "pCAMAL");
    strcpy(qa[0][1], "1.1");
    strncpy(qa[0][2], camal_date, MAX_STR_LENGTH);
    strncpy(qa[0][3], camal_time, MAX_STR_LENGTH);
    error = ex_put_qa(exoID, 1, qa);

    for (i = 0; i < 4; i++)
      delete [] qa[0][i];
  }

  return error;
}

int PCExodusFile::put_coor(int num_points_out, double *x_coor, 
                           double *y_coor, double *z_coor)
{
  int error = 1;
  
    // write coordinates
  if (num_points_out == numNodes) {
    error = ex_put_coord(exoID, x_coor, y_coor, z_coor);
  }

    // write coordinate names
  if (error == 0) {
    char *names[3] = {"x", "y", "z"};
    error = ex_put_coord_names(exoID, names);
  }

  return error;
}

int PCExodusFile::put_hex_blk(int blk_id, int num_nodes_elem,
                              int num_hexes, int *hexes)
{
  int num_attr = 1;
  int error = 1;
  char elem_type[8];
  switch (num_nodes_elem)
  {
    case 8:
        strcpy(elem_type, "HEX8");
        break;
    case 9:
    case 20:
    case 27:
        printf("WARNING: HEX%d element not supported in pCAMAL, using HEX8\n",
               num_nodes_elem);
        strcpy(elem_type, "HEX8");
        num_nodes_elem = 8;
        break;
    deault:
        printf("WARNING: Unknown element (%d-node hex), using HEX8\n",
                      num_nodes_elem);
        strcpy(elem_type, "HEX8");
        break;
  }
  
  
    // write block parameters
  if (num_hexes == numElems) {
    error = ex_put_elem_block(exoID, blk_id, elem_type, numElems, 
                              num_nodes_elem, num_attr);
  }
  
    // write element connectivity
  if (error == 0) {
    int i;
    for (i = 0; i < numElems * num_nodes_elem; i++)
      hexes[i] += 1; // increment node numbers by 1

    error = ex_put_elem_conn(exoID, blk_id, hexes);
  }
  
    // write element attributes
  if (error == 0) {
    double* attrib = new double[numElems * num_attr];
    if (attrib == NULL) {
      error = 1;
    }
    else {
      int i;
      for (i = 0; i < numElems * num_attr; i++) {
        attrib[i] = 1.0;
      }
      error = ex_put_elem_attr(exoID, blk_id, attrib);
      delete [] attrib;
    }
  }
  return error;
}

int  PCExodusFile::put_node_sets(int* node_id_array, int* ns_cnts_array, 
                                 int* ns_df_cnts_array, int* ns_ptrs, 
                                 int* ns_df_ptrs, int* ns_list,
                                 double* ns_df_list)
{
  if (exoID == 0 || numNodeSets == 0)
    return 1;
  
  return ex_put_concat_node_sets(exoID, node_id_array, ns_cnts_array, 
                                 ns_df_cnts_array, ns_ptrs, ns_df_ptrs,
                                 ns_list, ns_df_list);
}

void PCExodusFile::read_init()
{
    // open failed 
  if (exoID == 0)
    return;
  
    // read file init information  
  int error = ex_get_init(exoID, mTitle, &numDim, &numNodes, &numElems,
                          &numElemBlks, &numNodeSets, &numSideSets);

    // read element block IDs
  int eb_ids[numElemBlks];
  if (error == 0) {
    error = ex_get_elem_blk_ids(exoID, eb_ids);
  }

    // read "SweepID" property
  int surf_sweep1_ids[numElemBlks];
  memset(surf_sweep1_ids, 0, numElemBlks * sizeof(int));
  char *prop_name = "_CU_SweepID1";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, 
                              surf_sweep1_ids);
  }
  int surf_sweep2_ids[numElemBlks];
  memset(surf_sweep2_ids, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_SweepID2";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, 
                              surf_sweep2_ids);
  }

    // read "SurfaceType" property
  int surf_types1[numElemBlks];
  memset(surf_types1, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_SurfaceType1";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, surf_types1);
  }
  int surf_types2[numElemBlks];
  memset(surf_types2, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_SurfaceType2";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, surf_types2);
  }

    // read number of hexes that will be generated
  int num_hexes1[numElemBlks];
  memset(num_hexes1, 0, numElemBlks * sizeof(int));
  int num_hexes2[numElemBlks];
  memset(num_hexes2, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_NumberHexes1";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, num_hexes1);
  }
  if (error == 0) {
    prop_name = "_CU_NumberHexes2";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, num_hexes2);
  }
    // try old format if hexes not found
  else {
    prop_name = "_CU_NumberHexes";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, num_hexes1);    
  }

    // read user specified block id
  int block1_ids[numElemBlks];
  memset(block1_ids, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_UserBlockNo1";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, block1_ids);
  }
  int block2_ids[numElemBlks];
  memset(block2_ids, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_UserBlockNo2";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, block2_ids);
  }

    // read user specified element type
  int nodes_per_hex1[numElemBlks];
  memset(nodes_per_hex1, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_NodesPerHex1";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, nodes_per_hex1);
  }
  int nodes_per_hex2[numElemBlks];
  memset(nodes_per_hex2, 0, numElemBlks * sizeof(int));
  prop_name = "_CU_NodesPerHex2";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, nodes_per_hex2);
  }

    // generate sweep volumes for all SweepIDs
  std::map<int, PCSweepVolume*> sweep_map;
  if (error == 0) {
    delete_sweep_volumes();
    error = convert_sweep_data(eb_ids, block1_ids, surf_sweep1_ids, 
                               surf_types1, num_hexes1, nodes_per_hex1,
                               sweep_map);
  }
  if (error == 0) {
    error = convert_sweep_data(eb_ids, block2_ids, surf_sweep2_ids, 
                               surf_types2, num_hexes2, nodes_per_hex2,
                               sweep_map);
  }

    // update quad count (source, linking, target) for each sweep volume
  if (error == 0)
    error = update_quad_count();

  if (error != 0) {
    delete_sweep_volumes();
    numElemBlks = 0;
  }
}

int PCExodusFile::convert_sweep_data(int* eb_ids, int* block_ids, 
                                     int* surf_sweep_ids, int* surf_types, 
                                     int* num_hexes, int* nodes_per_hex,
                                     std::map<int, PCSweepVolume*>& sweep_map)
{
  PCSweepVolume* vol = NULL;
  int i;
  for (i = 0; i < numElemBlks; i++) {
    int sweep_id = surf_sweep_ids[i];

      // is this block a sweep block
    if (sweep_id <= 0)
      continue;
      
    if (sweep_map.empty() || sweep_map.find(sweep_id) == sweep_map.end()) {
      vol = new PCSweepVolume;
      sweepVols.push_back(vol);
      vol->put_sweep_id(sweep_id);
      if (block_ids[i] == 0)
        vol->put_elem_block_id(eb_ids[i]);  // no user defined block id
      else
        vol->put_elem_block_id(block_ids[i]);
      if (nodes_per_hex[i] == 0)
        vol->put_nodes_per_hex(8);  // default to 8-node hex
      else
        vol->put_nodes_per_hex(nodes_per_hex[i]);
      sweep_map[sweep_id] = vol;
    }
    else {
      vol = sweep_map[sweep_id];
    }

    switch (surf_types[i]) {
      case PCMLSweeper::SOURCE:
      case PCMLSweeper::TMP_SOURCE:
          vol->add_source_id(eb_ids[i]);
          vol->add_num_hexes(num_hexes[i]);
          break;
      case PCMLSweeper::LINKING:
      case PCMLSweeper::TMP_LINKING:
          vol->add_linking_id(eb_ids[i]);
          break;
      case PCMLSweeper::TARGET:
          vol->put_target_id(eb_ids[i]);
          break;
      default:
          printf("ERROR: Unknown surface type %d\n", surf_types[i]);
          break;
    }
  }
  return 0;
}

int PCExodusFile::update_quad_count()
{
    // nothing to update
  if (sweepVols.empty())
    return -1;
  
  std::set<int> blk_set;
  int error = 0;
  int i;
  for (i = 0; i < sweepVols.size() && error == 0; i++) {
    PCSweepVolume* vol = sweepVols[i];

    std::vector<int> blk_ids;
    vol->get_ordered_ids(blk_ids);
    int j;
    for (j = 0; j < blk_ids.size() && error == 0; j++) {
      char elem_type[MAX_STR_LENGTH+1];
      int num_elem_blk, num_node_elem, num_attr;
      int error = ex_get_elem_block(exoID, blk_ids[j], elem_type,
                                    &num_elem_blk, &num_node_elem, &num_attr);

      if (error == 0) {
        if (num_node_elem == 4)
          vol->put_num_surf_quads(num_elem_blk);
        else
          error = -1;
      }
    }
  }
  return error;
}

void PCExodusFile::delete_sweep_volumes()
{
  if (!sweepVols.empty()) {
    std::vector<PCSweepVolume*>::iterator it;
    for (it = sweepVols.begin(); it != sweepVols.end(); it++) {
      delete (*it);
    }
    sweepVols.clear();
  }
}

