#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"

#include "PCExodusFile.hpp"
#include "CMLSweeper.hpp"
#include "HexQuality.hpp"
#include "HexMeshQuality.hpp"

void WriteLocalExodusMesh( int num_points_out, int num_hexes, 
                           double* x_coor, double* y_coor, double* z_coor, 
			   int* connect, int sweep_id, char* fileout ) {
  // Write local Exodus II output file
  printf( "Writing an Exodus II mesh %d with %d points and %d hexes.\n",
	  sweep_id, num_points_out, num_hexes );

  char filename[strlen( fileout ) + 10];
  sprintf( filename, "%s.vol%03d.g", fileout, sweep_id );


  PCExodusFile exo_out( filename, pce::create );
  exo_out.put_param( num_points_out, num_hexes );
  exo_out.put_coor( num_points_out, x_coor, y_coor, z_coor );
  exo_out.put_hex_blk( num_hexes, connect );
}

int ReadSweepWriteSubdomains( PCExodusFile* pc_input, int vol_id, 
			      char* fileout, 
			      int& num_points_out, int& num_hexes,
                              double* q_mesh, int quality_measure, 
			      bool verbose ) {
  // Read sweep subdomain parameters
  int sweep_id, num_quads;
  pc_input->read_sweep_prop(vol_id, sweep_id, num_quads);
  if ( ! sweep_id )
    return 0;
  
  // Read coordinates
  int num_points;
  double* x_coor = NULL;
  double* y_coor = NULL;
  double* z_coor = NULL;
  int* node_ids  = NULL;
  pc_input->read_sweep_coord(vol_id, num_points, x_coor, y_coor, z_coor,
                             node_ids);
  if ( verbose )
    {
    printf( "\n			 ---Coordinates---\n" );
    printf( "  node	   x		y	     z\n" );
    int i;
    for (i = 0; i < num_points; i++) 
      {
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
  printf( "  number targets = %d\n", num_tgt_surf);
  printf( "	      total = %d\n", num_surfs);
  
  int num_surf_quads[num_surfs];
  pc_input->read_sweep_surf_size( vol_id, num_surfs, num_surf_quads );
  int num_tgt_quads = 0;
  for ( int i = 0; i < num_src_surf; ++i ) num_tgt_quads += num_surf_quads[i];
  printf( "\nNumber of quads/surface for subdomain %d\n", vol_id );
  printf( "  Surface	Quads\n" );
  for ( int i = 0; i < num_surfs; ++i ) 
    {
    printf( " %8d %8d\n", i+1, num_surf_quads[i]);
    }
  
  int connect[num_quads * 4];
  pc_input->read_sweep_conn( vol_id, num_points, num_quads, node_ids, connect );
  if ( verbose ) 
    {
    printf( "\n		  ---Connectivity---\n" );
    printf( "  Quad	n1	  n2	    n3	      n4\n" );
    int *c = connect;
    for ( int i = 0; i < num_quads; ++i ) 
      {
      printf( "%8d: %8d %8d %8d %8d\n", i+1, c[0], c[1], c[2], c[3]);
      c += 4;
      }
    }
  delete [] node_ids;

  // Setup CAMAL hex sweeper
  CMLSweeper sweeper;
  sweeper.set_boundary_mesh( num_points, x_coor, y_coor, z_coor,
			     num_quads, connect,
			     num_src_surf, num_surf_quads, num_tgt_quads );


  // Generate swept hex mesh
  sweeper.generate_mesh( num_points_out, num_hexes );
  delete [] x_coor;
  delete [] y_coor;
  delete [] z_coor;  

  // Retrieve mesh
  double x_coor_m[num_points_out];
  double y_coor_m[num_points_out];
  double z_coor_m[num_points_out];
  int connect_m[num_hexes * 8];
  sweeper.get_mesh( num_points_out, x_coor_m, y_coor_m, z_coor_m,
		    num_hexes, connect_m );

  // Get mesh quality
  HexMeshQuality hmq( x_coor_m, y_coor_m, z_coor_m, 
		      num_hexes, connect_m, quality_measure );
  q_mesh[0] = hmq.getMinQuality();
  q_mesh[1] = hmq.getMeanQuality();
  q_mesh[2] = hmq.getMaxQuality();
  q_mesh[3] = hmq.getStdvQuality();

  printf( "Mesh quality (min/mean/max/stdv) for subdomain %d: %G %G %G %G\n",
	  vol_id, q_mesh[0], q_mesh[1], q_mesh[2], q_mesh[3] );

  // Write mesh
  WriteLocalExodusMesh( num_points_out, num_hexes, 
			x_coor_m, y_coor_m, z_coor_m,
			connect_m, sweep_id, fileout );
  return 1;
}

void CalculateGlobalStats( int nproc, 
			   int* n_pts_g, int* n_pts_total, 
			   int* n_hex_g, int* n_hex_total,
                           double* q_mesh_g, double* q_mesh_total,
                           int verbose ) {
  int iproc;
  *n_pts_total = 0;
  *n_hex_total = 0;
  q_mesh_total[1] = 0.;
  
  if ( verbose )
    {
    printf( "\n# Local statistics:\n" );
    }
  for ( iproc = 0; iproc < nproc; ++ iproc )
    {
    if ( verbose )
      {
      printf( "   processor %d computed %d points and %d hexes\n", 
              iproc, n_pts_g[iproc], n_hex_g[iproc] );
      }
    *n_pts_total += n_pts_g[iproc];
    *n_hex_total += n_hex_g[iproc];
    q_mesh_total[1] += n_hex_g[iproc] * q_mesh_g[4 * iproc + 1];
    }
  q_mesh_total[1] /= *n_hex_total;
}

int BalanceLoad( int myrank, int nsub, int nproc, int* proc_assign, int verbose ) {

  if ( nproc < 1 ) 
    {
    if ( ! myrank ) 
      {
      printf( "## Error: not enough processors available.\n" );
      }
    return 1;
    }
  
  // FIXME: here, we will later distinguish between the 2 following subcases:
  //        1. nproc == nsub (unchanged)
  //        2. nproc > nsub (more work needed)
  if ( nproc >= nsub )
    {
    for ( int i = 0; i < nsub; ++ i ) proc_assign[i] = i;
    
    if ( ! myrank ) 
      {
      printf( "  load-balancing done (1 subdomain per processor).\n" );
      }
    return 0;
    }

  if ( nproc < nsub )
    {
    for ( int i = 0; i < nsub; ++ i ) proc_assign[i] = i % nproc;
    
    if ( ! myrank && verbose )
      {
      printf( "  load-balancing done:\n" );
      printf( "    subdomain  processor\n" );
      for ( int i = 0; i < nsub ; ++ i )
        {
        printf( "        %d         %d\n", i, proc_assign[i] ) ;
        }
      }
    else if ( ! myrank )
      {
      printf( "  load-balancing done (maximum of %d subdomains per processor).\n",
              (int) nsub / nproc );
      }

    return 0;
    }

  return 0;
}

int main(int argc, char **argv) {
  int myrank, nproc, nsub, nmesh;
  char filein[50];
  char fileout[50];
  bool verbose = false;

  if ( argc < 4 ) 
    {
    printf( "Usage: pcamal_proto <n_subdomains> <filein> <fileout>\n" );
    return 1;
    }

  nsub = atoi( argv[1] );
  strcpy( filein, argv[2] );
  strcpy( fileout, argv[3] );

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Comm_size( MPI_COMM_WORLD, &nproc );
  int n_pts_g[nproc], n_hex_g[nproc], proc_assign[nsub];
  double q_mesh_g[4*nproc];

  if ( ! myrank ) 
    {
    printf( "===============\n" );
    printf( "# Starting PCAMAL:\n", nsub );
    printf( "  number of subdomains: %d\n", nsub );
    printf( "  number of available processors: %d\n", nproc );
    }
  
  if ( BalanceLoad( myrank, nsub, nproc, proc_assign, 1 ) ) return 1;

  if ( ! myrank ) 
    {
    printf( "  input file name: %s\n", filein );
    }
  
  // Open input file
  PCExodusFile pc_input( filein, pce::read );
  int num_blks = pc_input.get_num_sweep_vols();

  if ( ! myrank ) printf( "\n# Input file read (%d subdomains).\n", num_blks );

  // Just for the sake of synchronizing printouts:
  MPI_Barrier( MPI_COMM_WORLD );

  int n_pts_l[1], n_hex_l[1];
  n_pts_l[0] = 0;
  n_hex_l[0] = 0;

  double q_mesh_l[4];
  q_mesh_l[0] = DBL_MAX; // min
  q_mesh_l[1] = 0.; // mean
  q_mesh_l[2] = 0.; // max
  q_mesh_l[3] = 0.; // stdv

  // Visit each subdomain
  for ( int vol_id = 0; vol_id < num_blks; ++ vol_id ) 
    {
    if ( myrank == proc_assign[vol_id] ) 
      {
      printf( "========\nProcess %d to handle subdomain %d.\n", myrank, vol_id );
      int num_points_out, num_hexes;
      double q_mesh[4];
      int sweepable = ReadSweepWriteSubdomains( &pc_input, vol_id, fileout, 
                                                num_points_out, num_hexes,
                                                q_mesh, PCAMAL_QUALITY_EDGE_RATIO,
						verbose );

      // Update local statistics
      q_mesh_l[1] = n_hex_l[0] * q_mesh_l[1] + num_hexes * q_mesh[1];
      n_pts_l[0] += num_points_out;
      n_hex_l[0] += num_hexes;
      q_mesh_l[1] /= n_hex_l[0];

      if ( q_mesh[0] < q_mesh_l[0] ) q_mesh_l[0] = q_mesh[0];

      if ( q_mesh[2] > q_mesh_l[2] ) q_mesh_l[2] = q_mesh[2];

      if ( ! sweepable ) 
        {
        printf( "## Subdomain %d was not to be swept.\n", vol_id );
        continue;
        }
      }
    }

  MPI_Gather( n_pts_l, 1, MPI_INTEGER, n_pts_g, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
  MPI_Gather( n_hex_l, 1, MPI_INTEGER, n_hex_g, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
  MPI_Gather( q_mesh_l, 4, MPI_DOUBLE, q_mesh_g, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  // Just for the sake of synchronizing printouts:
  MPI_Barrier( MPI_COMM_WORLD );

  if ( ! myrank ) 
    {
    int n_pts_total, n_hex_total;
    double q_mesh_total[4];
    CalculateGlobalStats( nproc, 
                          n_pts_g, &n_pts_total, 
                          n_hex_g, &n_hex_total, 
                          q_mesh_g, q_mesh_total, 1 );
    printf( "\n# Global statistics:\n" );
    printf( "  total number of points: %d\n", n_pts_total );
    printf( "  total number of hexes:  %d\n", n_hex_total );
    printf( "  Mesh quality (min/mean/max/stdv): %G %G %G %G\n",
	    q_mesh_total[0], q_mesh_total[1], q_mesh_total[2], q_mesh_total[3] );
    }
  
  MPI_Finalize();
  if ( ! myrank ) 
    {
    printf( "\n# Done. Exiting PCAMAL.\n" );
    printf( "===============\n\n" );
    }

  return 0;
}
