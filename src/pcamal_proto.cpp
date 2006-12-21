#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"

#include "PCExodusFile.hpp"
#include "CMLSweeper.hpp"

int ReadSweepWriteSubdomains( PCExodusFile* pc_input, int vol_id, 
			      char* fileout, bool verbose) {
  // Read sweep subdomain parameters
  int sweep_id, num_points, num_quads;
  pc_input->read_sweep_prop(vol_id, sweep_id, num_points, num_quads);
  if ( ! sweep_id )
    return 0;
  
  // Read coordinates
  double x_coor[num_points];
  double y_coor[num_points];
  double z_coor[num_points];
  int node_ids[num_points];
  pc_input->read_sweep_coord(vol_id, num_points, x_coor, y_coor, z_coor,
			    node_ids);
  if ( verbose ) {
    printf( "\n                  ---Coordinates---\n" );
    printf( "  node        x            y            z\n" );
    int i;
    for (i = 0; i < num_points; i++) {
      printf( "%8d %12.5e %12.5e %12.5e\n", 
	     i+1, x_coor[i], y_coor[i], z_coor[i]);
    }
  }
  
  // Read connectivity
  int num_src_surf, num_lnk_surf, num_tgt_surf;
  int num_surfs = pc_input->read_sweep_surf_prop(vol_id, num_src_surf, 
						num_lnk_surf, num_tgt_surf);
  printf( "\nSurface information for subdomain %d:\n", vol_id );
  printf( "  number sources = %d\n", num_src_surf);
  printf( "  number linking = %d\n", num_lnk_surf);
  printf( "  number target  = %d\n", num_tgt_surf);
  printf( "           total = %d\n", num_surfs);
  
  int num_surf_quads[num_surfs];
  pc_input->read_sweep_surf_size(vol_id, num_surfs, num_surf_quads);
  int num_tgt_quads = 0;
  for ( int i = 0; i < num_src_surf; ++i )
    num_tgt_quads += num_surf_quads[i];
  printf( "\nNumber of quads/surface for subdomain %d\n", vol_id );
  printf( "  Surface    Quads\n" );
  for ( int i = 0; i < num_surfs; ++i ) {
    printf( " %8d %8d\n", i+1, num_surf_quads[i]);
  }
  
  int connect[num_quads * 4];
  pc_input->read_sweep_conn(vol_id, num_quads, node_ids, connect);
  if ( verbose ) {
    printf( "\n           ---Connectivity---\n" );
    printf( "  Quad     n1        n2        n3        n4\n" );
    int *c = connect;
    for ( int i = 0; i < num_quads; ++i ) {
      printf( "%8d: %8d %8d %8d %8d\n", i+1, c[0], c[1], c[2], c[3]);
      c += 4;
    }
  }

  //CMLSweeper sweeper;
  //sweeper.set_boundary_mesh(num_points, x_coor, y_coor, z_coor,
  //			    num_quads, connect,
  //			    num_src_surf, num_surf_quads, num_tgt_quads);


  int num_points_out, num_hexes;

  // Write Exodus II output file
  char filename[strlen( fileout ) + 10];
  sprintf( filename, "%s.vol%03d.g", fileout, sweep_id );

  PCExodusFile exo_out(filename, pce::create);
  exo_out.put_param(num_points_out, num_hexes);
  exo_out.put_coor(num_points_out, x_coor, y_coor, z_coor);
  exo_out.put_hex_blk(num_hexes, connect);

  return 1;
}

void CalculateGlobalStats( int nsub, int* nglobal, double* qglobal, int* nmesh, double* qmesh ){
  int isub;
  *nmesh = 0;
  *qmesh = 0.;
  
//   for ( isub = 0; isub < nsub; ++ isub ) {
//     *nmesh += nglobal[isub];
//     *qmesh += qglobal[isub] * nglobal[isub];
//   }
//   *qmesh /= *nmesh;
}

int main(int argc, char **argv) {
  int myrank, nprocs, nsub, isub, nmesh;
  int nlocal[1], *nglobal;
  double qlocal[1], *qglobal, qmesh;
  char filein[50];
  char fileout[50];
  bool verbose = false;

  if ( argc < 4 ) {
    printf( "Usage: pcamal_proto <n_subdomains> <filein> <fileout>\n" );
    return 1;
  }

  nsub = atoi( argv[1] );
  strcpy( filein, argv[2] );
  strcpy( fileout, argv[3] );
  nglobal = (int*) malloc( nsub * sizeof( int ) );
  qglobal = (double*) malloc( nsub * sizeof( double ) );

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  if ( ! myrank ) {
    printf( "Number of subdomains: %d\n", nsub );
    printf( "Number of available processes: %d\n", nprocs );
    if ( nsub > nprocs ) {
      printf( "## Error: not enough processes available.\n" );
      return 1;
    }
    printf( "Input file to be read by all processes: %s\n", filein );
  }

  // Open input file
  PCExodusFile pc_input(filein, pce::read);
  int num_blks = pc_input.get_num_sweep_vols();

  MPI_Barrier( MPI_COMM_WORLD );
  if ( ! myrank ) printf( "Input file read (%d subdomains).\n", num_blks );

  // Visit each subdomain
  for ( int vol_id = 0; vol_id < num_blks; vol_id++) {
    if ( myrank == vol_id ) {
      printf( "Process %d to handle subdomain %d.\n", myrank, vol_id );
      int sweepable = ReadSweepWriteSubdomains( &pc_input, vol_id, fileout, false );
      if ( ! sweepable ) {
	printf( "## Subdomain %d was not to be swept.\n", vol_id );
	continue;
      }
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

  MPI_Finalize();
  return 0;
}
