#include "math.h"
#include "HexQuality.hpp"
#include "HexMeshQuality.hpp"

HexMeshQuality::HexMeshQuality( double* x_coor, double* y_coor, double* z_coor,
				int num_hexes, int* connect, int qualityMeasure )
{
  if ( num_hexes <= 0 ) 
    {
    this->min = this->max = this->mean = this->stdv = 0.;
    return;
    }

  double S2 = 0.;
  double q_elem;
  double coordinates[8][3];
  int *c = connect;
  for ( int j = 0; j < 8; ++j ) 
    {
    coordinates[j][0] = x_coor[c[j]];
    coordinates[j][1] = y_coor[c[j]];
    coordinates[j][2] = z_coor[c[j]];
    }
  q_elem = HexQuality::Shape( coordinates );
  this->min = this->max = this->mean = S2 = q_elem;
  c += 8;

  for ( int i = 1; i < num_hexes; ++i ) 
    {
    for ( int j = 0; j < 8; ++j ) 
      {
      coordinates[j][0] = x_coor[c[j]];
      coordinates[j][1] = y_coor[c[j]];
      coordinates[j][2] = z_coor[c[j]];
      }
    q_elem = HexQuality::Shape( coordinates );
    this->mean += q_elem;
    S2 += q_elem * q_elem;
    if ( q_elem < this->min) this->min = q_elem;
    if ( q_elem > this->max) this->max = q_elem;
    c += 8;
    }

  this->mean /= num_hexes;
  S2 /= num_hexes;
  this->stdv = sqrt( S2 - this->mean * this->mean );
}
