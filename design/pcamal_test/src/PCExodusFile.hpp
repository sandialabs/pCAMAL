// file: PCExodusFile.hpp
// description: header file for pCAMAL Exodus file interface
// author: Michael Stephenson

#ifndef PC_EXODUS_FILE_HPP
#define PC_EXODUS_FILE_HPP

#include <map>
#include <vector>
#include "exodusII.h"

namespace pce
{
    enum FileOp{create, read, write};
};


class PCSweepVolume;

class PCExodusFile
{
public:
  PCExodusFile(char* const filename, pce::FileOp op);
  virtual ~PCExodusFile();
  
  int  get_num_sweep_vols();
  int  get_num_hexes(int num_blks, int* num_hexes);
  void read_sweep_prop(int vol_id, int &sweep_id, int &num_quads);
  void read_sweep_coord(int vol_id, int& num_points, 
                        double* &x_coor, double* &y_coor, 
                        double* &z_coor, int* &node_ids);
  void read_sweep_conn(int vol_id, int num_points, int num_quads,
                       int *node_ids, int *connect);
  int  read_sweep_surf_prop(int vol_id, int &num_src_surf, int &num_lnk_surf,
                            int &num_tgt_surf);
  void read_sweep_surf_size(int vol_id, int num_surfs, int *num_surf_quads);

  int  put_param(int num_points_out, int num_hexes);
  int  put_coor(int num_points_out, double *x_coor, 
                double *y_coor, double *z_coor);
  int  put_hex_blk(int num_hexes, int *hexes);

private:
  int exoID;
  int cpuWord;
  int fileWord;
  int numDim;
  int numNodes;
  int numElems;
  int numElemBlks;
  int numNodeSets;
  int numSideSets;
  float mVersion;
  char mTitle[MAX_LINE_LENGTH+1];
  std::vector<PCSweepVolume*> sweepVols;

  int  convert_sweep_data(int* eb_ids, int* surf_sweep_ids, 
                          int* surf_types, int* num_hexes,
                          std::map<int, PCSweepVolume*>& sweep_map);
  void delete_sweep_volumes();
  void read_init();
  int  update_quad_count();
};


#endif // PC_EXODUS_FILE_HPP
