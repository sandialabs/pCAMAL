#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

void ExecuteCAMAL( int isub, int myrank, int* nlocal, double* qlocal ) {
  printf( "CAMAL called on subdomain %d by process %d:\n", isub, myrank );
  nlocal[0] = rint( drand48() * 10000 );
  qlocal[0] = drand48();
  printf( "    %d elements generated.\n", nlocal[0] );
  printf( "    average quality: %g\n", qlocal[0] );
}

void CalculateGlobalStats( int nsub, int* nglobal, double* qglobal, int* nmesh, double* qmesh ){
  int isub;
  *nmesh = 0;
  *qmesh = 0.;
  
  for ( isub = 0; isub < nsub; ++ isub ) {
    *nmesh += nglobal[isub];
    *qmesh += qglobal[isub] * nglobal[isub];
  }
  *qmesh /= *nmesh;
}

void SaveLocalFiles( int isub, int myrank, int nlocal ) {
  printf( "Process %d ready to save meshed subdomain %d:\n", myrank, isub );
  printf( "    %d elements to be saved.\n", nlocal );
}

int main(int argc, char **argv) {
  int myrank, nprocs, nsub, isub, nmesh;
  int nlocal[1], *nglobal;
  double qlocal[1], *qglobal, qmesh;
  char filename[30];

  if ( argc < 3 ) {
    printf("Usage: pcamal_proto <n_subdomains> <filename>\n");
    return 1;
  }

  nsub = atoi( argv[1] );
  strcpy( filename, argv[2] );
  nglobal = (int*) malloc( nsub * sizeof( int ) );
  qglobal = (double*) malloc( nsub * sizeof( double ) );

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  srand48( time( 0 ) * myrank );
  if ( ! myrank ) {
    printf( "Number of subdomains: %d\n", nsub );
    printf( "Number of available processes: %d\n", nprocs );
    if ( nsub > nprocs ) {
      printf("## Error: not enough processes available.\n");
      return 1;
    }
    printf( "Input file to be read from process 0: %s\n", filename );
  }

  // Just for the sake of synchronizing printouts:
  MPI_Barrier( MPI_COMM_WORLD );
  if ( ! myrank ) printf( "\n" );

  for ( isub = 0; isub < nsub; ++ isub ) {
    if ( isub == myrank ) {
      ExecuteCAMAL( isub, myrank, nlocal, qlocal );
    }
  }

  MPI_Gather( nlocal, 1, MPI_INTEGER, nglobal, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
  MPI_Gather( qlocal, 1, MPI_DOUBLE_PRECISION, qglobal, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD );

  // Just for the sake of synchronizing printouts:
  MPI_Barrier( MPI_COMM_WORLD );
  if ( ! myrank ) printf( "\n" );

  if ( ! myrank ) {
    CalculateGlobalStats( nsub, nglobal, qglobal, &nmesh, &qmesh );
    printf( "Global statistics:\n" );
    printf( "    total number of mesh elements: %d\n", nmesh );
    printf( "    average quality: %g\n", qmesh );
  }

  // Just for the sake of synchronizing printouts:
  MPI_Barrier( MPI_COMM_WORLD );
  if ( ! myrank ) printf( "\n" );

  for ( isub = 0; isub < nsub; ++ isub ) {
    if ( isub == myrank ) {
      SaveLocalFiles( isub, myrank, nlocal[0] );
    }
  }

  // Just for the sake of synchronizing printouts:
  MPI_Barrier( MPI_COMM_WORLD );
  if ( ! myrank ) printf( "\n" );

  MPI_Finalize();
  return 0;
}
