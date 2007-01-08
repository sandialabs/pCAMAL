// file: PCSweepVolume.cpp
// description: implementation of pCAMAL sweep volume
// author: Michael Stephenson

#include "exodusII.h"
#include "PCSweepVolume.hpp"

int PCSweepVolume::get_ordered_ids(std::vector<int>& ids)
{
  ids = sourceIds;
  int i;
  for (i = 0; i < linkingIds.size(); i++)
    ids.push_back(linkingIds[i]);
  ids.push_back(targetId);
  
  return ids.size();
}

int PCSweepVolume::get_num_quads()
{
  int num_quads = 0;
  int i;
  for (i = 0; i < numSurfQuads.size(); i++)
    num_quads += numSurfQuads[i];

  return num_quads;
}


int PCSweepVolume::get_num_surfs(int &num_src_surf, int &num_lnk_surf, 
                                 int &num_tgt_surf)
{
  num_src_surf = sourceIds.size();
  num_lnk_surf = linkingIds.size();
  num_tgt_surf = 1;
  
  return num_src_surf + num_lnk_surf + num_tgt_surf;
}

void PCSweepVolume::get_num_surf_quads(int num_surfs, int *num_surf_quads)
{
    // copy std::vector to int array
  int i;
  if (num_surfs == numSurfQuads.size()) {
    for (i = 0; i < num_surfs; i++) {
      num_surf_quads[i] = numSurfQuads[i];
    }
  }
  else {
    memset(num_surf_quads, 0, num_surfs * sizeof(int));
  }
}
