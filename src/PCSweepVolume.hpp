// file: PCSweepVolume.hpp
// description: header file for pCAMAL sweep volume
// author: Michael Stephenson

#ifndef PC_SWEEP_VOLUME_HPP
#define PC_SWEEP_VOLUME_HPP

#include <vector>

class PCSweepVolume
{
public:
  PCSweepVolume() : sweepId(0), numHexes(0), numQuads(0), elemBlockId(0), 
                    targetId(0), nodesPerHex(0)
      {}
  virtual ~PCSweepVolume()
      {}

  int  get_elem_block_id() {return elemBlockId;}
  void put_elem_block_id(int blk_id) {elemBlockId = blk_id;}

  int  add_num_hexes(int num_hexes) {numHexes += num_hexes;}
  int  get_num_hexes() {return numHexes;}
  void put_num_hexes(int num_hexes) {numHexes = num_hexes;}
  int  get_num_quads();

  int  get_num_surfs(int &num_src_surf, int &num_lnk_surf, int &num_tgt_surf);
  void get_num_surf_quads(int num_surfs, int *num_surf_quads);
  
  int  get_num_surf_quads(int i)
      { if (i >= 0 && i < numSurfQuads.size()) return numSurfQuads[i];
        else return 0;}
  void put_num_surf_quads(int num_quads) {numSurfQuads.push_back(num_quads);}

  int  get_sweep_id() {return sweepId;}
  void put_sweep_id(int sweep_id) {sweepId = sweep_id;}

  void add_source_id(int sid) {sourceIds.push_back(sid);}
  void add_linking_id(int lid) {linkingIds.push_back(lid);}
  void put_target_id(int tid) {targetId = tid;}
  int  get_ordered_ids(std::vector<int>& ids);

  void put_nodes_per_hex(int n) {nodesPerHex = n;}
  int  get_nodes_per_hex() {return nodesPerHex;}

private:
  int sweepId;
  int numHexes;
  int numQuads;
  int elemBlockId;
  int targetId;
  int nodesPerHex;
  std::vector<int> sourceIds;
  std::vector<int> linkingIds;
  std::vector<int> numSurfQuads;
};

#endif // PC_SWEEP_VOLUME_HPP

