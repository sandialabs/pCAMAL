#ifndef PC_HEX_QUALITY_HPP
#define PC_HEX_QUALITY_HPP

class PCHexQuality
{
public:
  static double EdgeRatio( double coordinates[][3] );
  static double Shape( double coordinates[][3] );
 
private:
  PCHexQuality() {};
  ~PCHexQuality() {};
};

#endif // PC_HEX_QUALITY_HPP
