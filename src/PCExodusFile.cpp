// file: PCExodusFile.cpp
// description: implementation of pCAMAL Exodus file interface
// author: Michael Stephenson

#include <map>
#include <string.h>
#include "PCExodusFile.hpp"
#include "PCSweepVolume.hpp"

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
  
int PCExodusFile::get_num_elem_blks()
{
  if (exoID == 0)
    return 0;

  return numElemBlks;
}

int PCExodusFile::get_num_sweep_vols()
{
  if (exoID == 0)
    return 0;

  return sweepVols.size();
}

void PCExodusFile::read_sweep_prop(int i, int &sweep_id, 
                                   int &num_points, int &num_quads)
{
  if (exoID == 0 || sweepVols.empty() ||
      i < 0 || i >= sweepVols.size()) {
    sweep_id   = 0;
    num_points = 0;
    num_quads  = 0;
  }
  else {
    sweep_id   = sweepVols[i]->get_sweep_id();
    num_points = sweepVols[i]->get_num_nodes();
    num_quads  = sweepVols[i]->get_num_quads();
  }  
}

void PCExodusFile::read_sweep_coord(int i, int num_points, 
                                    double *x_coor, double *y_coor, 
                                    double *z_coor, int *node_ids)
{
  if (num_points > 0 &&
      x_coor != NULL && y_coor != NULL && z_coor != NULL &&
      node_ids != NULL) {
    int ns_id = sweepVols[i]->get_node_set_id();
    int error = ex_get_node_set(exoID, ns_id, node_ids);

      // retrieve all nodes
    double *xx = NULL;
    double *yy = NULL;
    double *zz = NULL;
    if (error == 0) {
      xx = new double[numNodes];
      yy = new double[numNodes];
      zz = new double[numNodes];
      error = ex_get_coord(exoID, xx, yy, zz);      
    }

      // copy only coordinates in node set
    if (error == 0) {
      int i;
      for (i = 0; i < num_points; i++) {
        int j = node_ids[i] - 1;
        x_coor[i] = xx[j];
        y_coor[i] = yy[j];
        z_coor[i] = zz[j];
//         printf("%5d %5d: %12.5e %12.5e %12.5e\n",
//                i+1, j+1, xx[j], yy[j], zz[j]);
      }
    }

      // clean up
    delete [] zz;
    delete [] yy;
    delete [] xx;
  }
}

void PCExodusFile::read_sweep_conn(int i, int num_quads, int *node_ids,
                                   int *connect)
{
    // no Exodus file or no sweep volumes or index out of range
  if (exoID == 0 || sweepVols.empty() || i < 0 || i >= sweepVols.size())
    return;

    // read the connectivity
  if (num_quads > 0 && node_ids != NULL && connect != NULL) {
    int *tmp_conn = new int[num_quads * 4];
    int blk_id = sweepVols[i]->get_elem_block_id();
    int prev_blks = sweepVols[i]->get_elem_block_offset();    

    int error = ex_get_elem_conn(exoID, blk_id, tmp_conn);

      // sort elements by side set order
    int *c = connect;
    if (error == 0) {
      int nsrc, nlnk, ntgt;
      int num_surfs = sweepVols[i]->get_num_surfs(nsrc, nlnk, ntgt);
      int j;
      for (j = 0; j < num_surfs; j++) {
        int ssid = sweepVols[i]->get_side_set_id(j);
        int numq = sweepVols[i]->get_num_surf_quads(j);
        int *ss_list = new int[numq];
        int *ss_side = new int[numq];
        error = ex_get_side_set(exoID, ssid, ss_list, ss_side);

          // order tmp_conn in connect
        int k;
        for (k = 0; k < numq && error == 0; k++) {
          int index = ss_list[k];
          if (i == 0 && j == 1 && index > num_quads) 
            index -= num_quads;
          index = (index - prev_blks - 1) * 4;
          c[0] = tmp_conn[index++];
          c[1] = tmp_conn[index++];
          c[2] = tmp_conn[index++];
          c[3] = tmp_conn[index];
//           printf("%5d %5d: %5d %5d %5d %5d\n", 
//                  k, ss_list[k], c[0], c[1], c[2], c[3]);
          c += 4;
        }
        delete [] ss_side;
        delete [] ss_list;
      }
    }  
    delete [] tmp_conn;
    
      // map node numbers for node_ids to sequential ids
    if (error == 0) {
      int num_points = sweepVols[i]->get_num_nodes();
      std::map<int,int> node_map;
      int i;
      for (i = 0; i < num_points; i++) {
        node_map[node_ids[i]] = i;
      }
      
      for (i = 0; i < num_quads * 4; i++) {
        int id = connect[i];
        connect[i] = node_map[id];
//         printf("%5d: %5d -> %5d\n", i, id, connect[i]);
      }
    }
  }
}

int PCExodusFile::read_sweep_surf_prop(int i, int &num_src_surf,
                                        int &num_lnk_surf,
                                        int &num_tgt_surf)
{
  if (exoID == 0 || sweepVols.empty() ||
      i < 0 || i >= sweepVols.size()) {
    num_src_surf = 0;
    num_lnk_surf = 0;
    num_tgt_surf = 0;
  }
  else {
    int total = sweepVols[i]->get_num_surfs(num_src_surf, num_lnk_surf,
                                            num_tgt_surf);
  }
}

void PCExodusFile::read_sweep_surf_size(int i, int num_surfs,
                                        int *num_surf_quads)
{
  if (exoID == 0 || sweepVols.empty() ||
      i < 0 || i >= sweepVols.size())
    return;
  
  if (num_surfs > 0 && num_surf_quads != NULL) {
    sweepVols[i]->get_num_surf_quads(num_surfs, num_surf_quads);
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
    double *attrib = new double[numElems * num_attr];
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
  int *eb_ids = NULL;
  if (error == 0) {
    eb_ids = new int[numElemBlks];
    error = ex_get_elem_blk_ids(exoID, eb_ids);
  }

    // read node set IDs
  int *ns_ids = NULL;
  if (error == 0) {
    ns_ids = new int[numNodeSets];
    error = ex_get_node_set_ids(exoID, ns_ids);
  }
  
    // read node set IDs
  int *ss_ids = NULL;
  if (error == 0) {
    ss_ids = new int[numSideSets];
    error = ex_get_side_set_ids(exoID, ss_ids);
  }

    // read "SweepID" property
  int *blk_sweep_ids = NULL;
  char *prop_name = "SweepID";
  if (error == 0) {
    blk_sweep_ids = new int[numElemBlks];
    error = ex_get_prop_array(exoID, EX_ELEM_BLOCK, prop_name, blk_sweep_ids);
  }
  int *ns_sweep_ids = NULL;
  if (error == 0) {
    ns_sweep_ids = new int[numNodeSets];
    error = ex_get_prop_array(exoID, EX_NODE_SET, prop_name, ns_sweep_ids);
  }
  int *ss_sweep_ids = NULL;
  if (error == 0) {
    ss_sweep_ids = new int[numSideSets];
    error = ex_get_prop_array(exoID, EX_SIDE_SET, prop_name, ss_sweep_ids);
  }

    // read "SurfaceType" property
  int *surf_types = NULL;
  prop_name = "SurfaceType";
  if (error == 0) {
    surf_types = new int[numSideSets];
    error = ex_get_prop_array(exoID, EX_SIDE_SET, prop_name, surf_types);
  }

    // generate sweep volumes for all SweepIDs
  if (error == 0) {
    delete_sweep_volumes();
    
    int block_offset = 0;
    int i;
    for (i = 0; i < numElemBlks; i++) {
      char elem_type[MAX_STR_LENGTH+1];
      int num_elem_blk, num_node_elem, num_attr;
      error = ex_get_elem_block(exoID, eb_ids[i], elem_type, &num_elem_blk,
                                &num_node_elem, &num_attr);

      if (error == 0) {
        int sweep_id = blk_sweep_ids[i];

          // is this block a sweep block
        if (sweep_id > 0 && num_node_elem == 4) {
          PCSweepVolume* vol = new PCSweepVolume;
          sweepVols.push_back(vol);
          
          vol->put_sweep_id(sweep_id);
          vol->put_elem_block_id(eb_ids[i]);
          vol->put_elem_block_offset(block_offset);
          vol->put_num_quads(num_elem_blk);
          

            // find associated node set (assumes only one node set/elem block)
          int j;
          int ns_id = 0;
          for (j = 0; j < numNodeSets && error == 0; j++) {
            int ns_sweep_id = ns_sweep_ids[j];
            if (ns_sweep_id == sweep_id) {
              ns_id = ns_ids[j];
              vol->put_node_set_id(ns_id);
              int num_nodes, num_dist;
              error = ex_get_node_set_param(exoID, ns_id, &num_nodes,
                                            &num_dist);
              if (error == 0)
                vol->put_num_nodes(num_nodes);
              break;
            }
          }  

            // find associated side sets
          for (j = 0; j < numSideSets && error == 0; j++) {
            int ss_sweep_id = ss_sweep_ids[j];
            if (ss_sweep_id == sweep_id) {
              int ss_id = ss_ids[j];
              vol->put_side_set_id(ss_id, surf_types[j]);
              int num_quads, num_dist;
              error = ex_get_side_set_param(exoID, ss_id, &num_quads,
                                            &num_dist);
              if (error == 0)
                vol->put_num_surf_quads(num_quads);
            }
          }

            // sort the sides sets by surface type
          if (error == 0)
            vol->sort_surf();
        }
      }
      block_offset += num_elem_blk;
    }
  }

    // clean up
  delete [] surf_types;
  delete [] ss_sweep_ids;
  delete [] ns_sweep_ids;
  delete [] blk_sweep_ids;
  delete [] ss_ids;
  delete [] ns_ids;
  delete [] eb_ids;
}

void PCExodusFile::delete_sweep_volumes()
{
  if (!sweepVols.empty()) {
    std::vector<PCSweepVolume*>::iterator it;
    for (it = sweepVols.begin(); it != sweepVols.end(); it++) {
      delete (*it);
    }
    sweepVols.erase(sweepVols.begin(), sweepVols.end());
  }
}

