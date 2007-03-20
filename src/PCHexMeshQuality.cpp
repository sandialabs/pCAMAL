#include "PCHexQuality.hpp"
#include "PCHexMeshQuality.hpp"

#include <cmath>
#include <string>

using namespace std;

typedef double (*QualityType)( double[][3] );

const char* HexQualityNames[] =
{ 
  "Edge Ratio",
  "Shape"
};

PCHexMeshQuality::PCHexMeshQuality( double* x_coor, double* y_coor, double* z_coor,
				    int num_hexes, int* connect, 
				    int qualityIndex, string& qualityName )
{
  if ( num_hexes <= 0 ) 
    {
    this->min = this->max = this->mean = this->mom2 = 0.;
    return;
    }
  
  qualityName = HexQualityNames[qualityIndex];

  QualityType hexQuality;
  switch ( qualityIndex )
    {
    case PCAMAL_QUALITY_EDGE_RATIO:
      hexQuality = PCHexQuality::EdgeRatio;
      break;
    case PCAMAL_QUALITY_SHAPE:
      hexQuality = PCHexQuality::Shape;
      break;
    default:
      hexQuality = PCHexQuality::Shape;
    }

  double q_elem;
  double coordinates[8][3];
  int *c = connect;
  for ( int j = 0; j < 8; ++j ) 
    {
    coordinates[j][0] = x_coor[c[j]];
    coordinates[j][1] = y_coor[c[j]];
    coordinates[j][2] = z_coor[c[j]];
    }
  q_elem = hexQuality( coordinates );
  this->min = this->max = this->mean = q_elem;
  this->mom2 = q_elem * q_elem;
  c += 8;

  for ( int i = 1; i < num_hexes; ++i ) 
    {
    for ( int j = 0; j < 8; ++j ) 
      {
      coordinates[j][0] = x_coor[c[j]];
      coordinates[j][1] = y_coor[c[j]];
      coordinates[j][2] = z_coor[c[j]];
      }
    q_elem = hexQuality( coordinates );
    this->mean += q_elem;
    this->mom2 += q_elem * q_elem;
    if ( q_elem < this->min) this->min = q_elem;
    if ( q_elem > this->max) this->max = q_elem;
    c += 8;
    }

  this->mean /= num_hexes;
  this->mom2 /= num_hexes;
}
