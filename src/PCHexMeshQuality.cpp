#include "PCHexMeshQuality.hpp"

#include "verdict.h"

#include <iostream>

using namespace std;

typedef double (*QualityType)( int, double[][3] );

const char* HexQualityNames[] =
{ 
  "EdgeRatio",
  "AspectRatio",
  "RadiusRatio",
  "AspectFrobenius",
  "MedAspectFrobenius",
  "MaxAspectFrobenius",
  "MinAngle",
  "CollapseRatio",
  "MaxAngle",
  "Condition",
  "ScaledJacobian",
  "Shear",
  "RelativeSizeSquared",
  "Shape",
  "ShapeAndSize",
  "Distortion",
  "MaxEdgeRatio",
  "Skew",
  "Taper",
  "Volume",
  "Stretch",
  "Diagonal",
  "Dimension",
  "Oddy",
  "ShearAndSize",
  "Jacobian",
  "Warpage",
  "AspectGamma",
  "Area",
  "AspectBeta"
};

PCHexMeshQuality::PCHexMeshQuality( double* x_coor, double* y_coor, double* z_coor,
				    int num_hexes, int* connect, 
				    int& qualityIndex, string& qualityName )
{
  if ( num_hexes <= 0 ) 
    {
    this->min = this->max = this->mean = this->mom2 = 0.;
    return;
    }
  
  QualityType hexQuality;
  switch ( qualityIndex )
    {
    case PCAMAL_QUALITY_EDGE_RATIO:
      hexQuality = v_hex_edge_ratio;
      break;
    case PCAMAL_QUALITY_MED_ASPECT_FROBENIUS:
      hexQuality = v_hex_med_aspect_frobenius;
      break;
    case PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS:
      hexQuality = v_hex_max_aspect_frobenius;
      break;
    case PCAMAL_QUALITY_MAX_EDGE_RATIO:
      hexQuality = v_hex_max_edge_ratio;
      break;
    case PCAMAL_QUALITY_SKEW:
      hexQuality = v_hex_skew;
      break;
    case PCAMAL_QUALITY_TAPER:
      hexQuality = v_hex_taper;
      break;
    case PCAMAL_QUALITY_VOLUME:
      hexQuality = v_hex_volume;
      break;
    case PCAMAL_QUALITY_STRETCH:
      hexQuality = v_hex_stretch;
      break;
    case PCAMAL_QUALITY_DIAGONAL:
      hexQuality = v_hex_diagonal;
      break;
    case PCAMAL_QUALITY_DIMENSION:
      hexQuality = v_hex_dimension;
      break;
    case PCAMAL_QUALITY_ODDY:
      hexQuality = v_hex_oddy;
      break;
    case PCAMAL_QUALITY_CONDITION:
      hexQuality = v_hex_condition;
      break;
    case PCAMAL_QUALITY_JACOBIAN:
      hexQuality = v_hex_jacobian;
      break;
    case PCAMAL_QUALITY_SCALED_JACOBIAN:
      hexQuality = v_hex_scaled_jacobian;
      break;
    case PCAMAL_QUALITY_SHEAR:
      hexQuality = v_hex_shear;
      break;
    case PCAMAL_QUALITY_SHAPE:
      hexQuality = v_hex_shape;
      break;
    case PCAMAL_QUALITY_RELATIVE_SIZE_SQUARED:
      hexQuality = v_hex_relative_size_squared;
      break;
    case PCAMAL_QUALITY_SHAPE_AND_SIZE:
      hexQuality = v_hex_shape_and_size;
      break;
    case PCAMAL_QUALITY_SHEAR_AND_SIZE:
      hexQuality = v_hex_shear_and_size;
      break;
    case PCAMAL_QUALITY_DISTORTION:
      hexQuality = v_hex_distortion;
      break;
    default:
      cout << "Incorrect quality index ( "
           << qualityIndex 
           << " ), using max_aspect_frobenius instead.\n";
      qualityIndex = PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS;
      hexQuality = v_hex_max_aspect_frobenius;
      break;
    }

  qualityName = HexQualityNames[qualityIndex];

  double q_elem;
  double coordinates[8][3];
  int *c = connect;
  for ( int j = 0; j < 8; ++j ) 
    {
    coordinates[j][0] = x_coor[c[j]];
    coordinates[j][1] = y_coor[c[j]];
    coordinates[j][2] = z_coor[c[j]];
    }
  q_elem = hexQuality( 8, coordinates );
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
    q_elem = hexQuality( 8, coordinates );
    cout << q_elem << endl;
    this->mean += q_elem;
    this->mom2 += q_elem * q_elem;
    if ( q_elem < this->min) this->min = q_elem;
    if ( q_elem > this->max) this->max = q_elem;
    c += 8;
    }

  this->mean /= num_hexes;
  this->mom2 /= num_hexes;
}
