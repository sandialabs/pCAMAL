// file: PCExodusFile.cpp
// description: implementation of pCAMAL Exodus file interface
// author: Michael Stephenson

#include <set>
#include <string.h>
#include "PCExodusFile.hpp"
#include "PCSweepVolume.hpp"
#include "PCMLSweeper.hpp"

PCExodusFile::PCExodusFile(char* const filename, pce::FileOp op)
        : exoID(0), cpuWord(8), fileWord(0), numDim(0), numNodes(0),
          numElems(0), numElemBlks(0), numNodeSets(0), numSideSets(0),
          mVersion(0.0)
{
  switch(op) {
      // open exising file for read only
    case pce::read:
        exoID = ex_open(filename, EX_READ, &cpuWord, &fileWord,
                        &mVersion);
        read_init();
        break;

          // create a new file
    case pce::create:
        fileWord = 8;
        exoID = ex_create(filename, EX_CLOBBER, &cpuWord, &fileWord);
        break;

          // open existing file for read/write
    case pce::write:
        exoID = ex_open(filename, EX_WRITE, &cpuWord, &fileWord,
                        &mVersion);
        read_init();
        break;
    default:
        break;
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

void PCExodusFile::read_sweep_prop(int vol_id, int &sweep_id, int &num_quads)
{
  if (exoID == 0 || sweepVols.empty() ||
      vol_id < 0 || vol_id >= sweepVols.size()) {
    sweep_id   = 0;
    num_quads  = 0;
  }
  else {
    sweep_id   = sweepVols[vol_id]->get_sweep_id();
    num_quads  = sweepVols[vol_id]->get_num_quads();
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
  double xx[numNodes];
  double yy[numNodes];
  double zz[numNodes];
  int error = ex_get_coord(exoID, xx, yy, zz);      

    // copy only coordinates in element blocks
  int mark[numNodes];
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

int PCExodusFile::put_param(int num_points_out, int num_hexes)
{
  if (exoID == 0 || num_points_out == 0 || num_hexes == 0)
    return 1;

    // write file parameters
  char title[MAX_LINE_LENGTH+1];
  numDim      = 3;
  numNodes    = num_points_out;
  numElems    = num_hexes;
  numElemBlks = 1;
  numNodeSets = 0;
  numSideSets = 0;
  int error = ex_put_init(exoID, title, numDim, numNodes, numElems,
                          numElemBlks, numNodeSets, numSideSets);

    // write qa record
  if (error == 0) {
    char *qa[1][4];
    int i;
    for (i = 0; i < 4; i++)
      qa[0][i] = new char[MAX_STR_LENGTH+1];
    
    time_t result = time(NULL);
    struct tm *timeptr = localtime(&result);
    char camal_date[MAX_STR_LENGTH+1];
    strftime(camal_date, MAX_STR_LENGTH, "%D", timeptr);
    char camal_time[MAX_STR_LENGTH+1];
    strftime(camal_time, MAX_STR_LENGTH, "%H:%M:%S", timeptr);
    
      // write pCamal QA record
    sprintf(qa[0][0], "pCAMAL output: %s: %s", camal_date, camal_time);
    strcpy(qa[0][1], "1.0a");
    strcpy(qa[0][2], camal_date);
    strcpy(qa[0][3], camal_time);
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

int PCExodusFile::put_hex_blk(int num_hexes, int *hexes)
{
  int blk_id = 100;
  int num_nodes_elem = 8;
  int num_attr = 1;
  int error = 1;
  
    // write block parameters
  if (num_hexes == numElems) {
    error = ex_put_elem_block(exoID, blk_id, "HEX8", numElems, 
                              num_nodes_elem, num_attr);
  }
  
    // write element connectivity
  if (error == 0) {
    int i;
    for (i = 0; i < numElems * 8; i++)
      hexes[i] += 1; // increment node numbers by 1

    error = ex_put_elem_conn(exoID, blk_id, hexes);
  }
  
    // write element attributes
  if (error == 0) {
    double attrib[numElems * num_attr];
    int i;
    for (i = 0; i < numElems * num_attr; i++) {
      attrib[i] = 1.0;
    }
    error = ex_put_elem_attr(exoID, blk_id, attrib);
  }
  return error;
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
  char *prop_name = "_CU_SweepID1";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, 
                              surf_sweep1_ids);
  }
  int surf_sweep2_ids[numElemBlks];
  prop_name = "_CU_SweepID2";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, 
                              surf_sweep2_ids);
  }

    // read "SurfaceType" property
  int surf_types1[numElemBlks];
  prop_name = "_CU_SurfaceType1";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, surf_types1);
  }
  int surf_types2[numElemBlks];
  prop_name = "_CU_SurfaceType2";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, surf_types2);
  }

    // read number of hexes that will be generated
  int num_hexes[numElemBlks];
  prop_name = "_CU_NumberHexes";
  if (error == 0) {
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, num_hexes);
  }

    // generate sweep volumes for all SweepIDs
  std::map<int, PCSweepVolume*> sweep_map;
  if (error == 0) {
    delete_sweep_volumes();
    error = convert_sweep_data(eb_ids, surf_sweep1_ids, surf_types1, 
                               sweep_map);
  }
  if (error == 0) {
    error = convert_sweep_data(eb_ids, surf_sweep2_ids, surf_types2, 
                               sweep_map);
  }

    // update hex count for each sweep volume
  if (error == 0)
    error = update_hex_count(num_hexes, surf_sweep1_ids, surf_types1,
                             surf_sweep2_ids, surf_types2);

    // update quad count (source, linking, target) for each sweep volume
  if (error == 0)
    error = update_quad_count();
}

int PCExodusFile::convert_sweep_data(int* eb_ids, 
                                     int* surf_sweep_ids, int* surf_types,
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
      vol->put_elem_block_id(eb_ids[i]);
      sweep_map[sweep_id] = vol;
    }
    else {
      vol = sweep_map[sweep_id];
    }

    switch (surf_types[i]) {
      case PCMLSweeper::SOURCE:
      case PCMLSweeper::TMP_SOURCE:
          vol->add_source_id(eb_ids[i]);
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

int PCExodusFile::update_hex_count(int* num_hexes, int* surf_sweep_ids1,
                                   int* surf_types1, int* surf_sweep_ids2,
                                   int* surf_types2)
{
    // nothing to update
  if (sweepVols.empty())
    return -1;

    // initialize array for summation
  int num_vols = sweepVols.size();
  int total_hexes[num_vols];
  int error = 0;
  int i;
  for (i = 0; i < num_vols; i++) {
    total_hexes[i] = 0;
    sweepVols[i]->put_num_hexes(0);
  }

    // sum hexes from each source surface
  for (i = 0; i < numElemBlks && error == 0; i++) {
    int type1 = surf_types1[i];
    int type2 = surf_types2[i];
    if (type1 == PCMLSweeper::SOURCE ||
        type1 == PCMLSweeper::TMP_SOURCE) {
      if (type2 == PCMLSweeper::SOURCE ||
          type2 == PCMLSweeper::TMP_SOURCE) {
        error = -1;
      }
      else {
        int  id = surf_sweep_ids1[i] - 1;
        total_hexes[id] += num_hexes[i];
      }
    }
    if (type2 == PCMLSweeper::SOURCE ||
        type2 == PCMLSweeper::TMP_SOURCE) {
      if (type1 == PCMLSweeper::SOURCE ||
          type1 == PCMLSweeper::TMP_SOURCE) {
        error = -1;
      }
      else {
        int  id = surf_sweep_ids2[i] - 1;
        total_hexes[id] += num_hexes[i];
      }
    }
  }

    // update number of hexes in each volume
  for (i = 0; i < num_vols && error == 0; i++) {
    sweepVols[i]->put_num_hexes(total_hexes[i]);
  }
  
  return error;
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

