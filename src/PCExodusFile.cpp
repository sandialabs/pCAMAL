// file: PCExodusFile.cpp
// description: implementation of pCAMAL Exodus file interface
// author: Michael Stephenson

#include <set>
#include <string.h>
#include "PCExodusFile.hpp"
#include "PCSweepVolume.hpp"
#include "PCMLSweeper.hpp"

PCExodusFile::PCExodusFile(char* const file_name, pce::FileOp op)
        : zeroBased(true), exoID(0), cpuWord(8), fileWord(0), numDim(0), 
          numNodes(0), numElems(0), numElemBlks(0), numNodeSets(0), 
          numSideSets(0), numOutputBlocks(0), mVersion(0.0)
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

  memset(num_hexes, 0, num_blks * sizeof(int));

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

  num_blks      = numOutputBlocks;
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

bool PCExodusFile::get_node_sets(int& new_node_sets, int* ns_id_array, 
                                 int* ns_cnts_array, int* ns_df_cnts_array, 
                                 int* ns_ptrs, int* ns_df_ptrs,
                                 int* &ns_list, double* &ns_df_list)
{
  bool debug_flag = false;
  
    // fail if no file opened or no node sets
  if (exoID == 0 || numNodeSets == 0)
    return false;

  int ns_len = 0;
  int ns_df_len = 0;
  float fdum;
  char cdum;

    // get size of node set list
  int error = ex_inquire(exoID, EX_INQ_NS_NODE_LEN, &ns_len, &fdum, &cdum);

    // get size of distribution factors list
  if (error == 0) {
    error = ex_inquire(exoID, EX_INQ_NS_DF_LEN, &ns_df_len, &fdum, &cdum);
  }
  
    // allocate memory for node set and distribution factor lists
  if (error == 0 && ns_len > 0 && ns_df_len > 0) {
    ns_list = new int[ns_len];
    ns_df_list = new double[ns_df_len];

      // return if memory allocation failed
    if (ns_list == NULL || ns_df_list == NULL) {
      delete [] ns_list;
      delete [] ns_df_list;
      ns_list = NULL;
      ns_df_list = NULL;
      return false;
    }
    
      // get all node sets (concatenated)
    error = ex_get_concat_node_sets(exoID, ns_id_array, ns_cnts_array,
                                    ns_df_cnts_array, ns_ptrs, ns_df_ptrs,
                                    ns_list, ns_df_list);
  }

  if (debug_flag && error == 0) {
    printf("Before conversion\n");
    print_concat_node_sets(ns_id_array, ns_cnts_array, ns_df_cnts_array,
                           ns_ptrs, ns_df_ptrs, ns_list, ns_df_list);
  }

  return error == 0 ? true : false;
}

void PCExodusFile::print_concat_node_sets(int* ns_id_array, int* ns_cnts_array,
                                          int* ns_df_cnts_array, int* ns_ptrs,
                                          int* ns_df_ptrs, int* ns_list, 
                                          double* ns_df_list)
{
  bool debug_flag = false;
  
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

    // compute size of lists
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

    // additional output to compare coordinate positions
  if (debug_flag) {
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


bool PCExodusFile::get_side_sets(int vol_id, int& new_side_sets, int* num_el,
                                 int* &ss_conn, int* ss_id_array, 
                                 int* ss_cnts_array, int* ss_df_cnts_array, 
                                 int* ss_ptrs, int* ss_df_ptrs, int* &ss_list, 
                                 int* &ss_side_list, double* &ss_df_list)
{
  bool debug_flag = false;
  
    // fail if no file opened or no side sets
  if (exoID == 0 || numSideSets == 0)
    return false;

  int ss_len = 0;
  int ss_df_len = 0;
  float fdum;
  char cdum;

    // get size of side set (and side set side) list
  int error = ex_inquire(exoID, EX_INQ_SS_ELEM_LEN, &ss_len, &fdum, &cdum);

      // get size of distribution factors list
  if (error == 0) {
    error = ex_inquire(exoID, EX_INQ_SS_DF_LEN, &ss_df_len, &fdum, &cdum);
  }

      // allocate memory for side set list and side set side list
  if (error == 0 && ss_len > 0 && ss_df_len > 0) {
    ss_list = new int[ss_len];
    ss_side_list = new int[ss_len];
    ss_df_list = new double[ss_df_len];

      // return if memory allocation failed
    if (ss_list == NULL || ss_side_list == NULL || ss_df_list == NULL) {
      delete [] ss_list;
      delete [] ss_side_list;
      delete [] ss_df_list;
      ss_list = ss_side_list = NULL;
      ss_df_list = NULL;
      return false;
    }
    
      // get all side sets (concatenated)
    error = ex_get_concat_side_sets(exoID, ss_id_array, ss_cnts_array,
                                    ss_df_cnts_array, ss_ptrs, ss_df_ptrs,
                                    ss_list, ss_side_list, ss_df_list);
  }

  if (debug_flag && error == 0) {
    printf("Before conversion\n");
    print_concat_side_sets(ss_id_array, ss_cnts_array, ss_df_cnts_array,
                           ss_ptrs, ss_df_ptrs, ss_list, ss_side_list,
                           ss_df_list);
  }

    // get face connectivity for the side sets
  if (error == 0) {
    return convert_side_sets_to_quad_conn(vol_id, ss_cnts_array, ss_ptrs,
                                          ss_len, ss_list, num_el, ss_conn);
  }
  return false;
}

bool PCExodusFile::convert_side_sets_to_quad_conn(int vol_id, 
                                                  int* ss_cnts, int* ss_ptrs, 
                                                  int ss_len, int* ss_list, 
                                                  int* num_el, int* &ss_conn)
{
  int debug_flag = false;
  
    // fail in no open file or no element blocks
  if (exoID == 0 || numElemBlks == 0)
    return false;
  
    // select only those elements in this block
  int blk_ids[numElemBlks];
  int error = ex_get_elem_blk_ids(exoID, blk_ids);
  
    // get number of elements in each element block
  int elem_index[numElemBlks+1];
  memset(elem_index, 0, numElemBlks+1 * sizeof(int));
  if (error == 0) {
    int i;
    for (i = 0; i < numElemBlks && error == 0; i++) {
      char elem_type[MAX_STR_LENGTH+1];
      int nodes_per_elem;
      int num_attr;
      error = ex_get_elem_block(exoID, blk_ids[i], elem_type, &elem_index[i+1],
                                &nodes_per_elem, &num_attr);
    }
  }
  
    // generate running total (or index of where element block starts)
  if (error == 0) {
    int i, j;
    elem_index[0] = 1;
    for (i = 0; i < numElemBlks; i++) {
      elem_index[i+1] += elem_index[i];
    }
    
    if (debug_flag) {
      printf("element index/count\n");
      for (i = 0; i < numElemBlks; i++) {
        printf(" %5d", elem_index[i]);
      }
      printf("\n");
    }
    
      // count the number of side set elements in each element block
    int  elem_cnt[numElemBlks];
    memset(elem_cnt, 0, numElemBlks * sizeof(int));
    for (i = 0; i < ss_len; i++) {
      int elem_id = ss_list[i];
      int blk_id = find_blk(elem_id, elem_index);
      ++elem_cnt[blk_id];
    }

    if (debug_flag) {
      for (i = 0; i < numElemBlks; i++) {
        printf(" %5d", elem_cnt[i]);
      }
      printf("\n");
    }

      // use a set to identify only the element blocks 
      // with side set elements in the current sweep volume
    PCSweepVolume* blk_ptr = sweepVols[vol_id];
    std::vector<int> ordered_ids;
    blk_ptr->get_ordered_ids(ordered_ids);
    std::set<int> surf_ids_set(ordered_ids.begin(), ordered_ids.end());
    for (i = 0; i < numElemBlks; i++) {
      if (surf_ids_set.find(blk_ids[i]) == surf_ids_set.end())
        elem_cnt[i] = 0;
    }
    
    if (debug_flag) {
      for (i = 0; i < numElemBlks; i++) {
        printf(" %5d", elem_cnt[i]);
      }
      printf("\n");
    }
  
      // get the element connectivity associated with the side sets
    error = get_local_side_set_conn(ss_cnts, ss_ptrs, ss_len, ss_list, 
                                    elem_index, elem_cnt,
                                    blk_ids, num_el, ss_conn);
  }

  return error == 0 ? true : false;
}



void PCExodusFile::print_concat_side_sets(int* ss_id_array, int* ss_cnts_array,
                                          int* ss_df_cnts_array, int* ss_ptrs,
                                          int* ss_df_ptrs, int* ss_list, 
                                          int* ss_side_list,
                                          double* ss_df_list)
{
    // header info
  int i, j, k;
  for (i = 0; i < numSideSets; i++) {
    printf("side set %d\n", i+1);
    printf("\tside_set id     = %d\n", ss_id_array[i]);
    printf("\tnum sides       = %d\n", ss_cnts_array[i]);
    printf("\tnum dist_fact   = %d\n", ss_df_cnts_array[i]);
    printf("\tindex side_set  = %d\n", ss_ptrs[i]);
    printf("\tindex dist_fact = %d\n", ss_df_ptrs[i]);
  }

    // compute size of lists
  int ss_len = 0;
  int ss_df_len = 0;
  for (i = 0; i < numSideSets; i++) {
    ss_len += ss_cnts_array[i];
    ss_df_len += ss_df_cnts_array[i];
  }

    // side set
  if (ss_len > 0) {
    for (i = 0; i <numSideSets; i++) {
      printf("\nside_ss%d\n", i+1);
      int begin = ss_ptrs[i];
      int end   = begin + ss_cnts_array[i];

        // elements
      printf("element list\n");
      for (j = begin; j < end; j += 10) {
        printf("%4d:", j+1);
        for (k = 0; k < 10; k++) {
          if (j + k < end)
            printf(" %4d", ss_list[j+k]);
          else
            break;
        }
        printf("\n");
      }

        // sides
      printf("side list\n");
      for (j = begin; j < end; j += 10) {
        printf("%4d:", j+1);
        for (k = 0; k < 10; k++) {
          if (j + k < end)
            printf(" %4d", ss_side_list[j+k]);
          else
            break;
        }
        printf("\n");
      }
    }
  }

    // distribution factors
  if (ss_df_len > 0) {
    for (i = 0; i <numSideSets; i++) {
      printf("\ndist_fact_ss%d\n", i+1);
      int begin = ss_df_ptrs[i];
      int end   = begin + ss_df_cnts_array[i];
      for (j = begin; j < end; j += 5) {
        printf("%d:", j+1);
        for (k = 0; k < 5; k++) {
          if (j + k < end)
            printf(" %f", ss_df_list[j+k]);
          else
            break;
        }
        printf("\n");
      }
    }
  }
}

void PCExodusFile::print_side_sets_conn(int* num_el, int* ss_conn)
{
    // print the connectivity of all elements in the side set
  int num = 1;
  int* c = ss_conn;
  int begin = 0;
  int i;
  for (i = 0; i < numSideSets; i++) {
    int end = begin + num_el[i];
    int j;
    for (j = begin; j < end; j++) {
      printf("%5d: %5d %5d %5d %5d\n", num, c[0], c[1], c[2], c[3]);
      ++num;
      c += 4;
    }
  } 
}

int  PCExodusFile::get_local_side_set_conn(int* ss_cnts, int* ss_ptrs,
                                           int ss_len, int* ss_list, 
                                           int* elem_index, int* elem_cnt,
                                           int* blk_ids, int* num_el, 
                                           int* &ss_conn)
{
  bool debug_flag = false;
  int error = 0;
  int i, j;
  
    // get element connectivity for side sets
    // ignore element blocks that contain no side set elements
  int* blk_conn[numElemBlks];
  for (i = 0; i < numElemBlks; i++) {
    if (elem_cnt[i] == 0) {
      blk_conn[i] = NULL;
    }
    else {
      int* connect = new int[elem_cnt[i] * 4];
      error = ex_get_elem_conn(exoID, blk_ids[i], connect);
      blk_conn[i] = connect;
    }
  }

    // remove elements not in this volume using a set
  std::set<int> side_set_elems;
  for (i = 0; i < numElemBlks; i++) {
    if (elem_cnt[i] == 0) 
      continue;
    for (j = elem_index[i]; j < elem_index[i+1]; j++) {
      side_set_elems.insert(j);
    }
  }

    // mark elements not in this sweep volume with -1
  for (i = 0; i < ss_len; i++) {
    if (side_set_elems.find(ss_list[i]) == side_set_elems.end())
      ss_list[i] = -1;
  }
  if (debug_flag) {
    printf("\nSide set list before conversion\n");
    int begin = 0;
    for (i = 0; i < ss_len; i += 10) {
      int end = begin + 10;
      end = end < ss_len ? end : ss_len;
      printf("%5d:", i+1);
      for (j = begin; j < end; j++) {
        printf(" %5d", ss_list[j]);
      }
      printf("\n");
      begin = end;
    }
  }
  
  
    // copy side set connectivity to single output array
  if (error == 0) {
      // count number of element in side set and allocate memory
    int num_elems = 0;
    for (i = 0; i < ss_len; i++) {
      if (ss_len > 0)
        ++num_elems;
    }
    ss_conn = new int[num_elems * 4];
    if (ss_conn == NULL) // no memory
      return false;
    
    int kk = 1;
    int nn = 0;
    for (i = 0; i < numSideSets; i++) {
      num_el[i] = 0;
      int begin = ss_ptrs[i];
      int end   = begin + ss_cnts[i];
      for (j = begin; j < end; j++) {
        if (ss_list[j] > 0) {
          int elem_id = ss_list[j];
          int blk_id = find_blk(elem_id, elem_index);
          elem_id -= elem_index[blk_id];
          int* conn = blk_conn[blk_id];
          int* c = &conn[elem_id * 4];
          int k;
          for (k = 0; k < 4; k++) {
            ss_conn[nn++] = c[k];
          }
          ss_list[j] = kk++;
          ++num_el[i];
        }
      }
    }
  }

  if (debug_flag) {
    printf("\nSide set list after conversion\n");
    int begin = 0;
    for (i = 0; i < ss_len; i++) {
      if (ss_list[i] > 0) {
        int elem_id = ss_list[i] - 1;
        int* c = &ss_conn[elem_id * 4];
        printf("%5d - %5d: %5d %5d %5d %5d\n", i+1, ss_list[i],
               c[0], c[1], c[2], c[3]);
      }
      else {
        printf("%5d - %5d\n", i+1, ss_list[i]);
      }
    }
  }
  
    // delete connectivity memory
  for (i = 0; i < numElemBlks; i++) {
    if (blk_conn[i] != NULL)
      delete [] blk_conn[i];
  }

  return error;
}

int PCExodusFile::find_blk(int elem_id, int* elem_index)
{
    // use bisection search of ordered elem_index list
  int low  = 0;
  int high = numElemBlks;
  bool ascnd = elem_index[high] > elem_index[low];
  while (high - low > 1) {
    int mid = (high + low) >> 1;
    if (elem_id >= elem_index[mid] == ascnd)
      low = mid;
    else
      high = mid;
  }
  return low;
}

void PCExodusFile::read_sweep_coord(int vol_id, int& num_points, 
                                    double* &x_coor, double* &y_coor, 
                                    double* &z_coor, int* &node_ids)
{
  num_points = 0;
  if (vol_id < 0 || vol_id > sweepVols.size()) {
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
    
    int i, j;
    for (i = 0; i < blk_ids.size() && error == 0; i++) {
      int num_quads = vol->get_num_surf_quads(i);
      int conn[num_quads * 4];
      
      error = ex_get_elem_conn(exoID, blk_ids[i], conn);
      if (error == 0) {
        for (j = 0; j < num_quads * 4; j++) 
          mark[conn[j] - 1] = 1;
      }
    }
  }
    
    // count number of marked nodes
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
    increment_hexes(num_nodes_elem, hexes);
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

void PCExodusFile::increment_hexes(int num_nodes_elem, int* hexes)
{
    // convert zero-based to one-based connectivity ids
  if (zeroBased) {
    int i;
    for (i = 0; i < numElems * num_nodes_elem; i++)
      hexes[i] += 1; // increment node numbers by 1
  }
  zeroBased = false;
}

int PCExodusFile::put_node_sets(int num_points, int* node_ids, 
                                int* ns_id_array, int* ns_cnts_array, 
                                int* ns_df_cnts_array, int* ns_ptrs, 
                                int* ns_df_ptrs, int* ns_list,
                                double* ns_df_list)
{
  bool debug_flag = false;
  
    // fail if no file opened or no node sets
  if (exoID == 0 || numNodeSets == 0)
    return 1;
  
    // select only those nodes in this block (node_ids)
  int i, j;
  std::map<int,int> node_map;
  for (i = 0; i < num_points; i++) {
    node_map[node_ids[i]] = i + 1;
  }

    // compute size of lists
  int ns_len = 0;
  int ns_df_len = 0;
  for (i = 0; i < numNodeSets; i++) {
    ns_len += ns_cnts_array[i];
    ns_df_len += ns_df_cnts_array[i];
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

    // compute index and count arrays
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

    if (new_ns_len == ns_len) {
      printf("No change to side set\n");
    }
  }

  return ex_put_concat_node_sets(exoID, ns_id_array, ns_cnts_array, 
                                 ns_df_cnts_array, ns_ptrs, ns_df_ptrs,
                                 ns_list, ns_df_list);
}

int  PCExodusFile::put_side_sets(int num_points, int* node_ids,
                                 int* num_el, int* ss_conn, 
                                 int num_hexes, int* hexes,
                                 int* ss_id_array, int* ss_cnts_array, 
                                 int* ss_df_cnts_array, int* ss_ptrs, 
                                 int* ss_df_ptrs, int* ss_list, 
                                 int* ss_side_list, double* ss_df_list)
{
  bool debug_flag = false;
  
    // fail if no file opened or no side sets
  if (exoID == 0 || numSideSets == 0)  
    return 1;

  if (debug_flag) {
    printf("\nSide set connectivity before conversion\n");
    print_side_sets_conn(num_el, ss_conn);
  }

    // map new ids to old ids
  std::map<int,int> node_map;
  int i;
  for (i = 0; i < num_points; i++) {
    node_map[node_ids[i]] = i + 1;
  }
    
    // convert connectivity to new node ids
  int num_elems = 0;
  for (i = 0; i < numSideSets; i++) {
    num_elems += num_el[i];
  }
  int* c = ss_conn;
  int j;
  for (j = 0; j < num_elems * 4; j++) {
    if (node_map.find(c[j]) == node_map.end()) {
      c[j] = -1;
      printf("found bad node at index %d\n", j);
    }
    else {
      c[j] = node_map[c[j]];
    }
  }  

  if (debug_flag) {
    printf("\nSide set connectivity after conversion\n");
    print_side_sets_conn(num_el, ss_conn);
  }

    // increment zero-based connectivity
  const int num_nodes_elem = 8;
  increment_hexes(num_nodes_elem, hexes);

    // storage for hex element side set
  int hex_ss_ptrs[numSideSets];
  int hex_ss_df_ptrs[numSideSets];
  int hex_ss_df_cnts[numSideSets];
  int* hex_ss_cnts = num_el;
  int* hex_ss_list = new int[num_elems];
  int* hex_ss_side_list = new int[num_elems];
  double* hex_ss_df_list = new double[num_elems];
  if (hex_ss_list == NULL || hex_ss_side_list == NULL ||
      hex_ss_df_list == NULL) {
    delete [] hex_ss_df_list;
    delete [] hex_ss_side_list;
    delete [] hex_ss_list;
    printf("Error: failed to write side sets\n");
    return 1;
  }
  memset(hex_ss_list, 0, num_elems * sizeof(int));
  memset(hex_ss_side_list, 0, num_elems * sizeof(int));
  memset(hex_ss_df_list, 0, num_elems * sizeof(double));
  

  hex_ss_ptrs[0] = 0;
  hex_ss_df_ptrs[0] = 0;
  hex_ss_df_cnts[0] = hex_ss_cnts[0];
  for (i = 1; i < numSideSets; i++) {
    hex_ss_ptrs[i] += hex_ss_cnts[i-1];
    hex_ss_df_ptrs[i] += hex_ss_df_cnts[i-1];
    hex_ss_df_cnts[i] = hex_ss_cnts[i];
  }
  
    // generate list of hexes, sides and distribution factors for side sets
  const unsigned char code[8] =
      {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
  unsigned char marks[num_hexes];
  int first = 0;
  int ss_num;
  for (ss_num = 0; ss_num < numSideSets; ss_num++) {
    int last = first + num_el[ss_num] * 4;

      // put side set nodes in std::set and mark hexes
    std::set<int> side_set_nodes;
    for (i = first; i < last; i++) {
      assert(ss_conn[i] > 0);
      side_set_nodes.insert(ss_conn[i]);
    }

      // mark hexes with bit code indicating which node is in side set
    memset(marks, 0, num_hexes * sizeof(char));
    int* h = hexes;
    for (i = 0; i < num_hexes; i++) {
      for (j = 0; j < 8; j++) {
        assert(h[j] > 0);
        if (side_set_nodes.find(h[j]) != side_set_nodes.end()) {
          marks[i] |= code[j];
        }
      }
      h += 8;
    }

    if (debug_flag) {
      printf("\nHex codes\n");
      int begin = 0;
      for (i = 0; i < num_hexes; i += 10) {
        printf("%5d:", i+1);
        int end = begin + 10;
        end = end < num_hexes ? end : num_hexes;
        for (j = begin; j < end; j++) {
          printf(" %#5x", marks[j]);
        }
        printf("\n");
        begin = end;
      }
    }
    first = last;

      // copy element ids and side id
    int nn = hex_ss_ptrs[ss_num];
    for (i = 0; i < num_hexes; i++) {
      int side = 0;
      switch (marks[i]) {
        case 0x33: side = 1; break;
        case 0x66: side = 2; break;
        case 0xcc: side = 3; break;
        case 0x99: side = 4; break;
        case 0x0f: side = 5; break;
        case 0xf0: side = 6; break;
        default: break;
      }
      if (side > 0) {
        hex_ss_list[nn] = i+1;
        hex_ss_side_list[nn] = side;
        ++nn;
      }
    }
    assert(nn == hex_ss_cnts[ss_num]);

      // are distribution factors constant?
    bool constant_df = true;
    int begin = ss_df_ptrs[ss_num];
    int end   = begin + ss_df_cnts_array[ss_num];
    double df = ss_df_list[begin];
    for (i = begin + 1; i < end; i++) {
      if (df != ss_df_list[i]) {
        constant_df = false;
        break;
      }
    }

    if (!constant_df) {
      printf("Error: non-uniform distribution factors not supported yet!\n");
    }
    else {
      int begin = hex_ss_df_ptrs[ss_num];
      int end   = begin + hex_ss_df_cnts[ss_num];
      for (i = begin; i < end; i++) {
        hex_ss_df_list[i] = df;
      }
    }
  }

    // write the side sets
  int error = ex_put_concat_side_sets(exoID, ss_id_array, hex_ss_cnts,
                                      hex_ss_df_cnts, hex_ss_ptrs, 
                                      hex_ss_df_ptrs, hex_ss_list, 
                                      hex_ss_side_list, hex_ss_df_list);

    // delete memory
  delete [] hex_ss_df_list;
  delete [] hex_ss_side_list;
  delete [] hex_ss_list;

  return error;
}

void PCExodusFile::read_init()
{
  char* prop_name;
  
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
  if (error == 0) {
    prop_name = "_CU_SweepID1";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, 
                              surf_sweep1_ids);
  }
  int surf_sweep2_ids[numElemBlks];
  memset(surf_sweep2_ids, 0, numElemBlks * sizeof(int));
  if (error == 0) {
    prop_name = "_CU_SweepID2";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, 
                              surf_sweep2_ids);
  }

    // read "SurfaceType" property
  int surf_types1[numElemBlks];
  memset(surf_types1, 0, numElemBlks * sizeof(int));
  if (error == 0) {
    prop_name = "_CU_SurfaceType1";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, surf_types1);
  }
  int surf_types2[numElemBlks];
  memset(surf_types2, 0, numElemBlks * sizeof(int));
  if (error == 0) {
    prop_name = "_CU_SurfaceType2";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, surf_types2);
  }

    // read number of hexes that will be generated
  int num_hexes1[numElemBlks];
  memset(num_hexes1, 0, numElemBlks * sizeof(int));
  int num_hexes2[numElemBlks];
  memset(num_hexes2, 0, numElemBlks * sizeof(int));
  if (error == 0) {
    prop_name = "_CU_NumberHexes1";
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
  int block2_ids[numElemBlks];
  memset(block2_ids, 0, numElemBlks * sizeof(int));
  if (error == 0) {
    prop_name = "_CU_UserBlockNo1";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, block1_ids);
  }
  if (error == 0) {
    prop_name = "_CU_UserBlockNo2";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, block2_ids);
  }
    // nothing found (old format) so ignore
  else {
    error = 0;
  }

    // read user specified element type
  int nodes_per_hex1[numElemBlks];
  memset(nodes_per_hex1, 0, numElemBlks * sizeof(int));
  int nodes_per_hex2[numElemBlks];
  memset(nodes_per_hex2, 0, numElemBlks * sizeof(int));
  if (error == 0) {
    prop_name = "_CU_NodesPerHex1";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, nodes_per_hex1);
  }
  if (error == 0) {
    prop_name = "_CU_NodesPerHex2";
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, nodes_per_hex2);
  }
    // nothing found (old format) so ignore
  else {
    error = 0;
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

    // find other blocks (non-sweepable)
  int i;
  for (i = 0; i < numElemBlks; i++) {
    if (surf_sweep1_ids[i] == 0 && surf_sweep2_ids[i] == 0)
      ++numOutputBlocks;
  }
  numOutputBlocks += sweepVols.size();
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

