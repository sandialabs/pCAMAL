#ifndef PC_HEX_MESH_QUALITY_HPP
#define PC_HEX_MESH_QUALITY_HPP

#include <string>

class PCHexQuality;

#define PCAMAL_QUALITY_EDGE_RATIO 0
#define PCAMAL_QUALITY_SHAPE 1

class PCHexMeshQuality
{
public:
  PCHexMeshQuality( double* x_coor, double* y_coor, double* z_coor,
		    int num_hexes, int* connect, 
		    int qualityIndex, std::string& qualityName );

  virtual ~PCHexMeshQuality() {};

  double getMinQuality() const {return this->min;};
  double getMaxQuality() const {return this->max;};
  double getMeanQuality() const {return this->mean;};
  double getMom2Quality() const {return this->mom2;};

private:
  double min;
  double max;
  double mean;
  double mom2;
};

#endif // PC_HEX_MESH_QUALITY_HPP
