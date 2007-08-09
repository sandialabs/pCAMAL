#include "PCHexMeshQuality.hpp"

#include "verdict.h"

#include <cfloat>
#include <iostream>

using namespace std;

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

PCHexMeshQuality::PCHexMeshQuality( int qualID  )
{
  this->min = DBL_MAX;
  this->max = DBL_MIN;
  this->sum = this->sum2 = 0.;
  this->nHexes = this->minID = this->maxID = 0;
  this->qualID = qualID;
  
  switch ( qualID )
    {
    case PCAMAL_QUALITY_EDGE_RATIO:
      this->qualFct = v_hex_edge_ratio;
      break;
    case PCAMAL_QUALITY_MED_ASPECT_FROBENIUS:
      this->qualFct = v_hex_med_aspect_frobenius;
      break;
    case PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS:
      this->qualFct = v_hex_max_aspect_frobenius;
      break;
    case PCAMAL_QUALITY_MAX_EDGE_RATIO:
      this->qualFct = v_hex_max_edge_ratio;
      break;
    case PCAMAL_QUALITY_SKEW:
      this->qualFct = v_hex_skew;
      break;
    case PCAMAL_QUALITY_TAPER:
      this->qualFct = v_hex_taper;
      break;
    case PCAMAL_QUALITY_VOLUME:
      this->qualFct = v_hex_volume;
      break;
    case PCAMAL_QUALITY_STRETCH:
      this->qualFct = v_hex_stretch;
      break;
    case PCAMAL_QUALITY_DIAGONAL:
      this->qualFct = v_hex_diagonal;
      break;
    case PCAMAL_QUALITY_DIMENSION:
      this->qualFct = v_hex_dimension;
      break;
    case PCAMAL_QUALITY_ODDY:
      this->qualFct = v_hex_oddy;
      break;
    case PCAMAL_QUALITY_CONDITION:
      this->qualFct = v_hex_condition;
      break;
    case PCAMAL_QUALITY_JACOBIAN:
      this->qualFct = v_hex_jacobian;
      break;
    case PCAMAL_QUALITY_SCALED_JACOBIAN:
      this->qualFct = v_hex_scaled_jacobian;
      break;
    case PCAMAL_QUALITY_SHEAR:
      this->qualFct = v_hex_shear;
      break;
    case PCAMAL_QUALITY_SHAPE:
      this->qualFct = v_hex_shape;
      break;
    case PCAMAL_QUALITY_RELATIVE_SIZE_SQUARED:
      this->qualFct = v_hex_relative_size_squared;
      break;
    case PCAMAL_QUALITY_SHAPE_AND_SIZE:
      this->qualFct = v_hex_shape_and_size;
      break;
    case PCAMAL_QUALITY_SHEAR_AND_SIZE:
      this->qualFct = v_hex_shear_and_size;
      break;
    case PCAMAL_QUALITY_DISTORTION:
      this->qualFct = v_hex_distortion;
      break;
    default:
      cout << "Incorrect quality index ( "
           << this->qualID 
           << " ), using max_aspect_frobenius instead.\n";
      this->qualID = PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS;
      this->qualFct = v_hex_max_aspect_frobenius;
      break;
    }

  this->qualName = HexQualityNames[this->qualID];
}

int PCHexMeshQuality::Execute( int nHexes, double* x_coor, double* y_coor, double* z_coor, int* connect )
{
  if ( nHexes < 1 )
    {
    cout << "Incorrect number of elements ( "
         << nHexes 
         << " ), not calculating qualites.\n";
    return 1;
    }

  this->nHexes = nHexes;

  double q_elem;
  double coordinates[8][3];
  int *c = connect;
  for ( int j = 0; j < 8; ++j ) 
    {
    coordinates[j][0] = x_coor[c[j]];
    coordinates[j][1] = y_coor[c[j]];
    coordinates[j][2] = z_coor[c[j]];
    }

  q_elem = this->qualFct( 8, coordinates );
  this->min = this->max = this->sum = q_elem;
  this->sum2 = q_elem * q_elem;
  this->minID = this->maxID = 0;
  c += 8;

  for ( int i = 1; i < nHexes; ++i ) 
    {
    for ( int j = 0; j < 8; ++j ) 
      {
      coordinates[j][0] = x_coor[c[j]];
      coordinates[j][1] = y_coor[c[j]];
      coordinates[j][2] = z_coor[c[j]];
      }
    q_elem = this->qualFct( 8, coordinates );

    this->sum += q_elem;
    this->sum2 += q_elem * q_elem;
    if ( q_elem < this->min) 
      {
      this->min = q_elem;
      this->minID = i;
      }
    if ( q_elem > this->max) 
      {
      this->max = q_elem;
      this->maxID = i;
      }
    c += 8;
    }
  
  return 0;
}

std::string PCHexMeshQuality::qualIDToQualityName( int qualID )
{
  switch ( qualID )
    {
    case PCAMAL_QUALITY_EDGE_RATIO:
      return "EdgeRatio";
      break;
    case PCAMAL_QUALITY_MED_ASPECT_FROBENIUS:
      return "MedAspectFrobenius";
      break;
    case PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS:
      return "MaxAspectFrobenius";
      break;
    case PCAMAL_QUALITY_MAX_EDGE_RATIO:
      return "MaxEdgeRatio";
      break;
    case PCAMAL_QUALITY_SKEW:
      return "Skew";
      break;
    case PCAMAL_QUALITY_TAPER:
      return "Taper";
      break;
    case PCAMAL_QUALITY_VOLUME:
      return "Volume";
      break;
    case PCAMAL_QUALITY_STRETCH:
      return "Stretch";
      break;
    case PCAMAL_QUALITY_DIAGONAL:
      return "Diagonal";
      break;
    case PCAMAL_QUALITY_DIMENSION:
      return "Dimension";
      break;
    case PCAMAL_QUALITY_ODDY:
      return "Oddy";
      break;
    case PCAMAL_QUALITY_CONDITION:
      return "Condition";
      break;
    case PCAMAL_QUALITY_JACOBIAN:
      return "Jacobian";
      break;
    case PCAMAL_QUALITY_SCALED_JACOBIAN:
      return "ScaledJacobian";
      break;
    case PCAMAL_QUALITY_SHEAR:
      return "Shear";
      break;
    case PCAMAL_QUALITY_SHAPE:
      return "Shape";
      break;
    case PCAMAL_QUALITY_RELATIVE_SIZE_SQUARED:
      return "RelativeSizeSquared";
      break;
    case PCAMAL_QUALITY_SHAPE_AND_SIZE:
      return "ShapeAndSize";
      break;
    case PCAMAL_QUALITY_SHEAR_AND_SIZE:
      return "ShearAndSize";
      break;
    case PCAMAL_QUALITY_DISTORTION:
      return "Distortion";
      break;
    default:
      cout << "Incorrect quality index ( "
           << qualID 
           << " ), using MaxAspectFrobenius instead.\n";
      qualID = PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS;
      return "MaxAspectFrobenius";
      break;
    }
}
