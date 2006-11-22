// file: PCSweepVolume.hpp
// description: header file for pCAMAL sweep volume
// author: Michael Stephenson

#ifndef PC_SWEEP_VOLUME_HPP
#define PC_SWEEP_VOLUME_HPP

#include <vector>

class PCSweepVolume
{
public:
  PCSweepVolume() : sweepId(0),     numNodes(0),  numQuads(0), 
                    elemBlockId(0), nodeSetId(0), blockOffset(0)
      {}
  virtual ~PCSweepVolume()
      {}

  int  get_elem_block_id() {return elemBlockId;}
  void put_elem_block_id(int blk_id) {elemBlockId = blk_id;}
  
  int  get_node_set_id() {return nodeSetId;}
  void put_node_set_id(int set_id) {nodeSetId = set_id;}

  int  get_elem_block_offset() {return blockOffset;}
  void put_elem_block_offset(int offset) {blockOffset = offset;}

  int  get_side_set_id(int i)
      { if (i >= 0 && i < sideSetId.size()) return sideSetId[i];
        else return 0; }
  void put_side_set_id(int set_id, int type)
      {sideSetId.push_back(set_id); surfTypes.push_back(type);}

  int  get_num_nodes() {return numNodes;}
  void put_num_nodes(int num_nodes) {numNodes = num_nodes;}
  
  int  get_num_quads() {return numQuads;}
  void put_num_quads(int num_quads) {numQuads = num_quads;}

  int  get_num_surfs(int &num_src_surf, int &num_lnk_surf, int &num_tgt_surf);
  void get_num_surf_quads(int num_surfs, int *num_surf_quads);
  
  int  get_num_surf_quads(int i)
      { if (i >= 0 && i < numSurfQuads.size()) return numSurfQuads[i];
        else return 0;}
  void put_num_surf_quads(int num_quads) {numSurfQuads.push_back(num_quads);}

  int  get_sweep_id() {return sweepId;}
  void put_sweep_id(int sweep_id) {sweepId = sweep_id;}

  int  sort_surf();

private:
  int sweepId;
  int numNodes;
  int numQuads;
  int elemBlockId;
  int nodeSetId;
  int blockOffset;
  std::vector<int> sideSetId;
  std::vector<int> surfTypes;
  std::vector<int> numSurfQuads;
};

#endif // PC_SWEEP_VOLUME_HPP

