#ifndef HEX_MESH_QUALITY_HPP
#define HEX_MESH_QUALITY_HPP

class HexQuality;

#define PCAMAL_QUALITY_EDGE_RATIO 0
#define PCAMAL_QUALITY_SHAPE 1

class HexMeshQuality
{
public:
  HexMeshQuality( double* x_coor, double* y_coor, double* z_coor,
				int num_hexes, int* connect, int qualityMeasure );

  ~HexMeshQuality() {};

  double getMinQuality() const {return this->min;};
  double getMaxQuality() const {return this->max;};
  double getMeanQuality() const {return this->mean;};
  double getStdvQuality() const {return this->stdv;};

private:
  double min;
  double max;
  double mean;
  double stdv;
};

#endif // HEX_MESH_QUALITY_HPP
