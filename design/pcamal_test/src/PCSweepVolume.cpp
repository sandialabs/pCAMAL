// file: PCSweepVolume.cpp
// description: implementation of pCAMAL sweep volume
// author: Michael Stephenson

#include "exodusII.h"
#include "PCSweepVolume.hpp"

int PCSweepVolume::sort_surf()
{
    // check consistency of data
  if (sideSetId.size() == surfTypes.size() &&
      surfTypes.size() == numSurfQuads.size()) {

      // insertion sort since surfaces should be nearly in-order
    int i, j;
    for (j = 1; j < sideSetId.size(); j++) {
      int type = surfTypes[j];
      int ssid = sideSetId[j];
      int numQ = numSurfQuads[j];
      i = j - 1;
      while (i >= 0 && surfTypes[i] > type) {
        surfTypes[i+1] = surfTypes[i];
        sideSetId[i+1] = sideSetId[i];
        numSurfQuads[i+1] = numSurfQuads[i];
        --i;
      }
      surfTypes[i+1] = type;
      sideSetId[i+1] = ssid;
      numSurfQuads[i+1] = numQ;
    }
    return 0;
  }
  return -1;
}

int PCSweepVolume::get_num_surfs(int &num_src_surf, int &num_lnk_surf, 
                                 int &num_tgt_surf)
{
  num_src_surf = 0;
  num_lnk_surf = 0;
  num_tgt_surf = 0;
  
    // check consistency of data
  if (sideSetId.size() == surfTypes.size() &&
      surfTypes.size() == numSurfQuads.size()) {
    int i;
    for (i = 0; i < surfTypes.size(); i++) {
      switch(surfTypes[i]) {
        case 1:
            ++num_src_surf;
            break;
        case 2:
            ++num_lnk_surf;
            break;
        case 3:
            ++num_tgt_surf;
            break;
        default:
            break;
      }
    }
  }
  return num_src_surf + num_lnk_surf + num_tgt_surf;
}

void PCSweepVolume::get_num_surf_quads(int num_surfs, int *num_surf_quads)
{
    // check consistency of data
  if (sideSetId.size() == surfTypes.size() &&
      surfTypes.size() == numSurfQuads.size()) {

      // copy std::vector to int array
    if (num_surfs == numSurfQuads.size()) {
      int i;
      for (i = 0; i < num_surfs; i++) {
        num_surf_quads[i] = numSurfQuads[i];
      }
    }
  }
}
