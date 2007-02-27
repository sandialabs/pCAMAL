#include "math.h"
#include "stdio.h"
#include "HexQuality.hpp"

inline double length_squared( const double edge[3] )
{
  return edge[0] * edge[0] + edge[1] * edge[1] + edge[2] * edge[2];
}

inline void vectorDifference( const double coordinates[8][3], 
                              int a, int b, double vd[3] )
{
  for ( int i = 0; i < 3; ++ i )
    {
    vd[i] = coordinates[b][i] - coordinates[a][i];
    }
}

inline double det3x3( const double c1[3], 
                      const double c2[3], 
                      const double c3[3] )
{
  return c1[0]*c2[1]*c3[2] + c2[0]*c3[1]*c1[2] + c3[0]*c1[1]*c2[2] -
         c1[0]*c3[1]*c2[2] - c2[0]*c1[1]*c3[2] - c3[0]*c2[1]*c1[2];
}


inline void make_hex_edges( const double coordinates[8][3], 
                            double edges[12][3] )
{

  for ( int i = 0; i < 3; ++ i )
    {
    edges[0][i] = coordinates[1][i] - coordinates[0][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[1][i] = coordinates[2][i] - coordinates[1][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[2][i] = coordinates[3][i] - coordinates[2][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[3][i] = coordinates[0][i] - coordinates[3][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[4][i] = coordinates[5][i] - coordinates[4][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[5][i] = coordinates[6][i] - coordinates[5][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[6][i] = coordinates[7][i] - coordinates[6][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[7][i] = coordinates[4][i] - coordinates[7][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[8][i] = coordinates[4][i] - coordinates[0][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[9][i] = coordinates[5][i] - coordinates[1][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[10][i] = coordinates[6][i] - coordinates[2][i];
    }

  for ( int i = 0; i < 3; ++ i )
    {
    edges[11][i] = coordinates[7][i] - coordinates[3][i];
    }
}

double HexQuality::EdgeRatio( double coordinates[8][3] )
{
  double edges[12][3];
  make_hex_edges(coordinates, edges);
  
  double a2 = length_squared( edges[0] );
  double b2 = length_squared( edges[1] );
  double c2 = length_squared( edges[2] );
  double d2 = length_squared( edges[3] );
  double e2 = length_squared( edges[4] );
  double f2 = length_squared( edges[5] );
  double g2 = length_squared( edges[6] );
  double h2 = length_squared( edges[7] );
  double i2 = length_squared( edges[8] );
  double j2 = length_squared( edges[9] );
  double k2 = length_squared( edges[10] );
  double l2 = length_squared( edges[11] );

  double mab,mcd,mef,Mab,Mcd,Mef;
  double mgh,mij,mkl,Mgh,Mij,Mkl;

  if ( a2 < b2 )
    {
      mab = a2;
      Mab = b2;
    }
  else // b2 <= a2
    {
      mab = b2;
      Mab = a2;
    }
  if ( c2 < d2 )
    {
      mcd = c2;
      Mcd = d2;
    }
  else // d2 <= c2
    {
      mcd = d2;
      Mcd = c2;
    }
  if ( e2 < f2 )
    {
      mef = e2;
      Mef = f2;
    }
  else // f2 <= e2
    {
      mef = f2;
      Mef = e2;
    }
  if ( g2 < h2 )
    {
      mgh = g2;
      Mgh = h2;
    }
  else // h2 <= g2
    {
      mgh = h2;
      Mgh = g2;
    }
  if ( i2 < j2 )
    {
      mij = i2;
      Mij = j2;
    }
  else // j2 <= i2
    {
      mij = j2;
      Mij = i2;
    }
  if ( k2 < l2 )
    {
      mkl = k2;
      Mkl = l2;
    }
  else // l2 <= k2
    {
      mkl = l2;
      Mkl = k2;
    }

  double m2;
  m2 = mab < mcd ? mab : mcd;
  m2 = m2  < mef ? m2  : mef;
  m2 = m2  < mgh ? m2  : mgh;
  m2 = m2  < mij ? m2  : mij;
  m2 = m2  < mkl ? m2  : mkl;

  double M2;
  M2 = Mab > Mcd ? Mab : Mcd;
  M2 = M2  > Mef ? M2  : Mef;
  M2 = M2  > Mgh ? M2  : Mgh;
  M2 = M2  > Mij ? M2  : Mij;
  M2 = M2  > Mkl ? M2  : Mkl;
  m2 = m2  < mef ? m2  : mef;

  M2 = Mab > Mcd ? Mab : Mcd;
  M2 = M2  > Mef ? M2  : Mef;

  return sqrt( M2 / m2 );
}

double HexQuality::Shape( double coordinates[8][3] )
{
  double det, shape;
  double min_shape = 1.; 
  static const double two_thirds = 2. / 3.;

  double  xxi[3], xet[3], xze[3];

  // J(0,0,0):
  vectorDifference( coordinates, 1, 0, xxi );
  vectorDifference( coordinates, 3, 0, xet );
  vectorDifference( coordinates, 4, 0, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }
  

  // J(1,0,0):
  vectorDifference( coordinates, 2, 1, xxi );
  vectorDifference( coordinates, 0, 1, xet );
  vectorDifference( coordinates, 5, 1, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }
  

  // J(1,1,0):
  vectorDifference( coordinates, 3, 2, xxi );
  vectorDifference( coordinates, 1, 2, xet );
  vectorDifference( coordinates, 6, 2, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }
  

  // J(0,1,0):
  vectorDifference( coordinates, 0, 3, xxi );
  vectorDifference( coordinates, 2, 3, xet );
  vectorDifference( coordinates, 7, 3, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }
  

  // J(0,0,1):
  vectorDifference( coordinates, 7, 4, xxi );
  vectorDifference( coordinates, 5, 4, xet );
  vectorDifference( coordinates, 0, 4, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }
  

  // J(1,0,1):
  vectorDifference( coordinates, 4, 5, xxi );
  vectorDifference( coordinates, 6, 5, xet );
  vectorDifference( coordinates, 1, 5, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }
  

  // J(1,1,1):
  vectorDifference( coordinates, 5, 6, xxi );
  vectorDifference( coordinates, 7, 6, xet );
  vectorDifference( coordinates, 2, 6, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }
  

  // J(1,1,1):
  vectorDifference( coordinates, 6, 7, xxi );
  vectorDifference( coordinates, 4, 7, xet );
  vectorDifference( coordinates, 3, 7, xze );

  det = det3x3( xxi, xet, xze );
  shape = 3 * pow( det, two_thirds) / ( length_squared( xxi ) + length_squared( xet ) + length_squared( xze ) );
  
  if( shape < min_shape ) { min_shape = shape; }

  return min_shape;
}
