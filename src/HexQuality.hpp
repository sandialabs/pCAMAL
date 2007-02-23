#ifndef HEX_QUALITY_HPP
#define HEX_QUALITY_HPP

class HexQuality
{
public:
  static double EdgeRatio( double coordinates[][3] );

protected:
  HexQuality();
  ~HexQuality();

private:
  HexQuality( const HexQuality& ); // Not implemented.
  void operator = ( const HexQuality& ); // Not implemented.
};

#endif // HEX_QUALITY_HPP
