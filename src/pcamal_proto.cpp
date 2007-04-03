#include "mpi.h"

#include <cfloat>
#include <cmath>
#include <string>
#include <iostream>

#include "PCExodusFile.hpp"
#include "PCMLSweeper.hpp"
#include "PCHexMeshQuality.hpp"

using namespace std;

string qualityName; 

void WriteLocalExodusMesh( int num_points_out, int num_hexes, 
                           double* x_coor, double* y_coor, double* z_coor, 
			   int* connect, int sweep_id, char* fileout,
                           int verbose ) {
  // Write local Exodus II output file
  if ( verbose )
    {
    cout << "Writing Exodus II mesh " 
         << sweep_id 
         << " with " 
         << num_points_out 
         << " points and " 
         << num_hexes 
         << " hexes." << endl;
    }

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
                              double* q_mesh, int qualityIndex, 
			      int verbose ) {
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
  if ( verbose > 1 )
    {
    cout <<  "\n			 ---Coordinates---" << endl;
    cout <<  "  node	   x		y	     z" << endl;
    int i;
    for (i = 0; i < num_points; i++) 
      {
      cout << " "
	   << i+1 
	   << x_coor[i] 
	   << y_coor[i]
	   << z_coor[i]
	   << endl;
      }
    }
  
  // Read connectivity
  int num_src_surf, num_lnk_surf, num_tgt_surf;
  int num_surfs = pc_input->read_sweep_surf_prop(vol_id, num_src_surf, 
                                                 num_lnk_surf, num_tgt_surf);
  if ( verbose )
    {
    cout <<  endl 
         << "Surface information for subdomain " 
         << vol_id
         << ":"
         << endl
         <<  "  number sources = "
         << num_src_surf
         << endl
         <<  "  number linking = "
         << num_lnk_surf
         << endl
         <<  "  number targets = "
         << num_tgt_surf
         << endl
         <<  "	   total = " 
         << num_surfs
         << endl;
    }
  
  int num_surf_quads[num_surfs];
  pc_input->read_sweep_surf_size( vol_id, num_surfs, num_surf_quads );
  int num_tgt_quads = 0;
  for ( int i = 0; i < num_src_surf; ++i ) num_tgt_quads += num_surf_quads[i];
  if ( verbose )
    {
    cout <<  endl 
         << "Number of quads/surface for subdomain "
         << vol_id 
         << endl
         <<  "Surface  Quads" 
         << endl;
    for ( int i = 0; i < num_surfs; ++i ) 
      {
      cout << " "
	   << i+1
	   << ": "
	   << num_surf_quads[i]
	   << endl;
      }
    }
  
  int connect[num_quads * 4];
  pc_input->read_sweep_conn( vol_id, num_points, num_quads, node_ids, connect );
  if ( verbose > 1 ) 
    {
    cout <<  endl 
	 << "		  ---Connectivity---" 
	 << endl;
    cout <<  "  Quad	n1	  n2	    n3	      n4" 
	 << endl;
    int *c = connect;
    for ( int i = 0; i < num_quads; ++i ) 
      {
      cout << i+1
	   << ": "
	   << c[0]
	   << " "
	   << c[1]
	   << " "
	   << c[2]
	   << " "
	   << c[3]
	   << endl;
      c += 4;
      }
    }
  delete [] node_ids;

  // Setup CAMAL hex sweeper
  PCMLSweeper sweeper;
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
  PCHexMeshQuality hmq( x_coor_m, y_coor_m, z_coor_m, 
			num_hexes, connect_m, 
			qualityIndex, qualityName );
  q_mesh[0] = hmq.getMinQuality();
  q_mesh[1] = hmq.getMeanQuality();
  q_mesh[2] = hmq.getMaxQuality();
  q_mesh[3] = hmq.getMom2Quality();

    if ( verbose > 1 )
    {
    cout <<  "Mesh quality ("
         << qualityName.c_str()
         << ") of subdomain "
         << vol_id
         << ":"
         << " min= "
         << q_mesh[0]
         << " mean= "
         << q_mesh[1]
         << " max= "
         << q_mesh[2]
         << " stdv= "
         << sqrt( q_mesh[3] - q_mesh[1] * q_mesh[1] )
         << endl;
    }

  // Write mesh
  WriteLocalExodusMesh( num_points_out, num_hexes, 
			x_coor_m, y_coor_m, z_coor_m,
			connect_m, sweep_id, fileout,
                        verbose );
  return 1;
}

void CalculateGlobalStats( int nproc, 
			   int* n_pts_g, int& n_pts_total, 
			   int* n_hex_g, int& n_hex_total,
                           double* q_mesh_g, double* q_mesh_total,
                           int verbose ) {
  n_pts_total = 0;
  n_hex_total = 0;
  q_mesh_total[0] = DBL_MAX; // min
  q_mesh_total[1] = 0.; // mean
  q_mesh_total[2] = 0.; // max
  q_mesh_total[3] = 0.; // mom2
  
  if ( verbose )
    {
    cout << endl 
	 << "# Local statistics:" 
	 << endl;
    }
  int ix4  = 0;
  for ( int iproc = 0; iproc < nproc; ++ iproc, ix4 += 4 )
    {
    if ( verbose )
      {
      cout <<  "   processor "
	   << iproc
	   << " computed "
	   << n_pts_g[iproc]
	   << " points and "
	   << n_hex_g[iproc]
	   << " hexes"
	   << endl;
      }
    n_pts_total += n_pts_g[iproc];
    n_hex_total += n_hex_g[iproc];
    q_mesh_total[1] += n_hex_g[iproc] * q_mesh_g[ix4 + 1];
    q_mesh_total[3] += n_hex_g[iproc] * q_mesh_g[ix4 + 3];

    if ( q_mesh_g[ix4] < q_mesh_total[0] ) q_mesh_total[0] = q_mesh_g[ix4];
    
    if ( q_mesh_g[ix4 + 2] > q_mesh_total[2] ) q_mesh_total[2] = q_mesh_g[ix4 + 2];
    }
  q_mesh_total[1] /= n_hex_total;
  q_mesh_total[3] /= n_hex_total;
}

int BalanceLoad( int myrank, int nsub, int nproc, int* proc_assign, int verbose ) {

  if ( nproc < 1 ) 
    {
    if ( ! myrank ) 
      {
      cout <<  "## Error: not enough processors available." << endl;
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
      cout <<  "  load-balancing done (1 subdomain per processor)." << endl;
      }
    return 0;
    }

  if ( nproc < nsub )
    {
    for ( int i = 0; i < nsub; ++ i ) proc_assign[i] = i % nproc;
    
    if ( ! myrank && verbose )
      {
      cout <<  "  load-balancing done:" << endl;
      cout <<  "    subdomain  processor" << endl;
      for ( int i = 0; i < nsub ; ++ i )
        {
        cout <<  "        "
	     << i
	     << "         "
	     << proc_assign[i]
	     << endl;
        }
      }
    else if ( ! myrank )
      {
      cout <<  "  load-balancing done (maximum of "
	   << (int) nsub / nproc
	   << " subdomains per processor)."
	   << endl;
      }

    return 0;
    }

  return 0;
}

int main(int argc, char **argv) {
  int myrank, nproc, nsub, nmesh;
  char filein[50];
  char fileout[50];
  bool verbose = 0;

  if ( argc < 4 ) 
    {
    cout <<  "Usage: pcamal_proto <n_subdomains> <filein> <fileout>" << endl;
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
    cout <<  "===============" 
	 << endl
         <<  "# Starting PCAMAL:"
	 << endl
         <<  "  number of subdomains: "
	 << nsub
	 << endl
         <<  "  number of available processors: "
	 << nproc
	 << endl;
    }
  
  if ( BalanceLoad( myrank, nsub, nproc, proc_assign, 0 ) ) return 1;

  if ( ! myrank ) 
    {
    cout <<  "  input file name: "
	 << filein
	 << endl;
    }
  
  // Open input file
  PCExodusFile pc_input( filein, pce::read );
  int num_blks = pc_input.get_num_sweep_vols();

  if ( ! myrank ) 
    {
    cout <<  endl 
         << "# Input file read ("
         << num_blks
         << " subdomains)."
         << endl;
    }

  // Just for the sake of synchronizing printouts:
  MPI_Barrier( MPI_COMM_WORLD );

  int n_pts_l[1], n_hex_l[1];
  n_pts_l[0] = 0;
  n_hex_l[0] = 0;

  double q_mesh_l[4];
  q_mesh_l[0] = DBL_MAX; // min
  q_mesh_l[1] = 0.; // mean
  q_mesh_l[2] = 0.; // max
  q_mesh_l[3] = 0.; // mom2

  // Visit each subdomain
  for ( int vol_id = 0; vol_id < num_blks; ++ vol_id ) 
    {
    if ( myrank == proc_assign[vol_id] ) 
      {
      if ( verbose )
        {
        cout <<  "========"
             << endl
             << "Process "
             << myrank
             << " to handle subdomain "
             <<  vol_id
             << "."
             << endl;
        }
      int num_points_out, num_hexes;
      double q_mesh[4];
      int sweepable = ReadSweepWriteSubdomains( &pc_input, vol_id, fileout, 
                                                num_points_out, num_hexes,
                                                q_mesh, PCAMAL_QUALITY_SHAPE,
						verbose );

      // Update local statistics
      // mean & mom2
      q_mesh_l[1] = n_hex_l[0] * q_mesh_l[1] + num_hexes * q_mesh[1];
      q_mesh_l[3] = n_hex_l[0] * q_mesh_l[3] + num_hexes * q_mesh[3];
      n_pts_l[0] += num_points_out;
      n_hex_l[0] += num_hexes;
      q_mesh_l[1] /= n_hex_l[0];
      q_mesh_l[3] /= n_hex_l[0];

      // min & max
      if ( q_mesh[0] < q_mesh_l[0] ) q_mesh_l[0] = q_mesh[0];
      if ( q_mesh[2] > q_mesh_l[2] ) q_mesh_l[2] = q_mesh[2];

      if ( ! sweepable ) 
        {
        cout <<  "## Subdomain "
             << vol_id
             << "was not to be swept."
             << endl;
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
                          n_pts_g, n_pts_total, 
                          n_hex_g, n_hex_total, 
                          q_mesh_g, q_mesh_total, 
			  1 );
    cout <<  endl 
         << "# Global statistics:" 
         << endl
         <<  "  total number of points: "
         << n_pts_total
         << endl
         <<  "  total number of hexes: "
         << n_hex_total
         << endl
         <<  "  mesh quality ("
         << qualityName.c_str()
         << "):"
         << " min= "
         << q_mesh_total[0]
         << " mean= "
         << q_mesh_total[1]
         << " max= "
         << q_mesh_total[2]
         << " stdv= "
         << sqrt( q_mesh_total[3] - q_mesh_total[1] * q_mesh_total[1] )
         << endl;
    }
  
  MPI_Finalize();
  if ( ! myrank ) 
    {
    cout <<  endl 
	 << "# Done. Exiting PCAMAL." 
	 << endl
	 <<  "===============" 
	 << endl;
    }

  return 0;
}
