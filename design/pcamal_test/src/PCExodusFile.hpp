// file: PCExodusFile.hpp
// description: header file for pCAMAL Exodus file interface
// author: Michael Stephenson

#ifndef PC_EXODUS_FILE_HPP
#define PC_EXODUS_FILE_HPP

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
  
  int get_num_elem_blks();
  int get_num_sweep_vols();
  void read_sweep_prop(int i, int &sweep_id, 
                       int &num_points, int &num_quads);
  void read_sweep_coord(int i, int num_points, 
                        double *x_coor, double *y_coor, 
                        double *z_coor, int *node_ids);
  void read_sweep_conn(int i, int num_quads, int *node_ids,
                       int *connect);
  int  read_sweep_surf_prop(int i, int &num_src_surf, int &num_lnk_surf,
                            int &num_tgt_surf);
  void read_sweep_surf_size(int i, int num_surfs, int *num_surf_quads);

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

  void delete_sweep_volumes();
  void read_init();
};


#endif // PC_EXODUS_FILE_HPP
