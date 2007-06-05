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
  void get_param(int& num_blks, int& num_node_sets, int& num_side_sets);
  bool get_node_sets(int& new_node_sets, int* ns_id_array, int* ns_cnts_array,
                     int* ns_df_cnts_array, int* ns_ptrs, int* ns_df_ptrs,
                     int* &ns_list, double* &ns_df_list);
  bool get_side_sets(int vol_id, int& new_side_sets, int* num_el,
                     int* &ss_conn, int* ss_id_array, int* ss_cnts_array, 
                     int* ss_df_cnts_array, int* ss_ptrs, int* ss_df_ptrs, 
                     int* &ss_list, int* &ss_side_list, double* &ss_df_list);

  void read_sweep_prop(int vol_id, int &block_id, int &sweep_id, 
                       int &num_quads, int& nodes_per_hex);
  void read_sweep_coord(int vol_id, int& num_points, 
                        double* &x_coor, double* &y_coor, 
                        double* &z_coor, int* &node_ids);
  void read_sweep_conn(int vol_id, int num_points, int num_quads,
                       int *node_ids, int *connect);
  int  read_sweep_surf_prop(int vol_id, int &num_src_surf, int &num_lnk_surf,
                            int &num_tgt_surf);
  void read_sweep_surf_size(int vol_id, int num_surfs, int *num_surf_quads);

  int  put_param(int num_points_out, int num_hexes, 
                 int num_node_sets, int num_side_sets);
  int  put_coor(int num_points_out, double *x_coor, 
                double *y_coor, double *z_coor);
  int  put_hex_blk(int block_id, int num_nodes_elem, 
                   int num_hexes, int *hexes);
  int  put_node_sets(int num_points, int* node_ids, 
                     int* node_id_array, int* ns_cnts_array, 
                     int* ns_df_cnts_array, int* ns_ptrs, int* ns_df_ptrs, 
                     int* ns_list, double* ns_df_list);
  int  put_side_sets(int num_points, int* node_ids, int* elem_cnt, 
                     int* ss_conn, int num_hexes, int* hexes,
                     int* side_id_array, int* ss_cnts_array, 
                     int* ss_df_cnts_array, int* ss_ptrs, int* ss_df_ptrs, 
                     int* ss_list, int* ss_side_list, double* ss_df_list);

private:
  bool zeroBased;
  int exoID;
  int cpuWord;
  int fileWord;
  int numDim;
  int numNodes;
  int numElems;
  int numElemBlks;
  int numNodeSets;
  int numSideSets;
  int numOutputBlocks;
  float mVersion;
  char fileName[MAX_STR_LENGTH+1];
  char mTitle[MAX_LINE_LENGTH+1];
  std::vector<PCSweepVolume*> sweepVols;

  int  convert_sweep_data(int* eb_ids, int* block_ids, int* surf_sweep_ids, 
                          int* surf_types, int* num_hexes, int* nodes_per_hex,
                          std::map<int, PCSweepVolume*>& sweep_map);

  bool convert_side_sets_to_quad_conn(int vol_id, int* ss_cnts, int* ss_ptrs,
                                      int ss_len, int* ss_list, 
                                      int* num_el, int* &ss_conn);
  
  void delete_sweep_volumes();

  int  find_blk(int elem_id, int* elem_index);

  void increment_hexes(int num_nodes_elem, int* hexes);
  
  int  get_local_side_set_conn(int* ss_cnts, int* ss_ptrs,
                               int ss_len, int* ss_list,
                               int* elem_index, int* elem_cnt,
                               int* blk_ids, int* num_el, int* &ss_conn);
  
  void print_concat_node_sets(int* ns_id_array, int* ns_cnts_array, 
                              int* ns_df_cnts_array, int* ns_ptrs, 
                              int* ns_df_ptrs, int* ns_list, 
                              double* ns_df_list);
  void print_concat_side_sets(int* ss_id_array, int* ss_cnts_array,
                              int* ss_df_cnts_array, int* ss_ptrs,
                              int* ss_df_ptrs, int* ss_list, 
                              int* ss_side_list,
                              double* ss_df_list);
  void print_node_sets_coords(int* ns_id_array, int* ns_cnts_array,
                              int* ns_ptrs, int* ns_list,
                              double* xx, double* yy, double* zz);
  void print_side_sets_conn(int* elem_cnt, int* ss_conn);

  void read_init();
  int  update_quad_count();
};


#endif // PC_EXODUS_FILE_HPP
