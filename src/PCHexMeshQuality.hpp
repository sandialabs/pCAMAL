#ifndef PC_HEX_MESH_QUALITY_HPP
#define PC_HEX_MESH_QUALITY_HPP

#include <string>

#define PCAMAL_QUALITY_EDGE_RATIO 0
#define PCAMAL_QUALITY_ASPECT_RATIO 1
#define PCAMAL_QUALITY_RADIUS_RATIO 2
#define PCAMAL_QUALITY_ASPECT_FROBENIUS 3
#define PCAMAL_QUALITY_MED_ASPECT_FROBENIUS 4
#define PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS 5
#define PCAMAL_QUALITY_MIN_ANGLE 6
#define PCAMAL_QUALITY_COLLAPSE_RATIO 7
#define PCAMAL_QUALITY_MAX_ANGLE 8
#define PCAMAL_QUALITY_CONDITION 9
#define PCAMAL_QUALITY_SCALED_JACOBIAN 10
#define PCAMAL_QUALITY_SHEAR 11
#define PCAMAL_QUALITY_RELATIVE_SIZE_SQUARED 12
#define PCAMAL_QUALITY_SHAPE 13
#define PCAMAL_QUALITY_SHAPE_AND_SIZE 14
#define PCAMAL_QUALITY_DISTORTION 15
#define PCAMAL_QUALITY_MAX_EDGE_RATIO 16
#define PCAMAL_QUALITY_SKEW 17
#define PCAMAL_QUALITY_TAPER 18
#define PCAMAL_QUALITY_VOLUME 19
#define PCAMAL_QUALITY_STRETCH 20
#define PCAMAL_QUALITY_DIAGONAL 21
#define PCAMAL_QUALITY_DIMENSION 22
#define PCAMAL_QUALITY_ODDY 23
#define PCAMAL_QUALITY_SHEAR_AND_SIZE 24
#define PCAMAL_QUALITY_JACOBIAN 25
#define PCAMAL_QUALITY_WARPAGE 26
#define PCAMAL_QUALITY_ASPECT_GAMMA 27
#define PCAMAL_QUALITY_AREA 28
#define PCAMAL_QUALITY_ASPECT_BETA 29

typedef double (*QualityType)( int, double[][3] );

class PCHexMeshQuality
{
public:
  PCHexMeshQuality( int qualID );

  int Execute( int nHexes, double* x_coor, double* y_coor, double* z_coor, int* connect );

  virtual ~PCHexMeshQuality() {};

  double getMinQuality() const {return this->min;};
  double getMaxQuality() const {return this->max;};
  double getSumValuesQuality() const {return this->sum;};
  double getSumSquaresQuality() const {return this->sum2;};
  int getNumberHexes() const {return this->nHexes;};
  int getQualityID() const {return this->qualID;};
  int getMinQualityID() const {return this->minID;};
  int getMaxQualityID() const {return this->maxID;};
  std::string getQualityName() const {return this->qualName;};

  static std::string qualIDToQualityName( int qualID );

private:
  double min;
  double max;
  double sum;
  double sum2;

  int nHexes;
  int qualID;
  int minID;
  int maxID;

  QualityType qualFct;

  std::string qualName;
};

#endif // PC_HEX_MESH_QUALITY_HPP
