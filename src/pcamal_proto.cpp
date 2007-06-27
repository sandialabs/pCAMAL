#include "mpi.h"
#include "sys/time.h"

#include <cfloat>
#include <cmath>
#include <string>
#include <iostream>

#include "PCExodusFile.hpp"
#include "PCMLSweeper.hpp"
#include "PCHexMeshQuality.hpp"

using namespace std;

string qualityName; 

void ConvertToPatranOrder(int num_hexes, int* connect)
{
  int tmp;
  int *c = connect;
  int i;
  for (i = 0; i < num_hexes; i++) {
    tmp  = c[2];
    c[2] = c[5];
    c[5] = tmp;
    tmp  = c[3];
    c[3] = c[4];
    c[4] = tmp;
    c += 8;
  }
}

void WriteLocalExodusMesh( PCExodusFile* pc_input, int vol_id,
                           int num_points_in, int* node_ids,
                           int num_points_out, int num_hexes, 
                           int num_node_sets, int num_side_sets,
                           double* x_coor, double* y_coor, double* z_coor, 
                           int* connect, int sweep_id, int blk_id,
                           int num_nodes_elem, char* fileout, int verbose ) {
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
  exo_out.put_param( num_points_out, num_hexes, num_node_sets, num_side_sets );

    // read/modify/write node set if any from input file
  if (num_node_sets > 0)
  {
    int new_node_sets = 0;
    int ns_len = 0;
    int ns_df_len = 0;
    int ns_ids[num_node_sets];
    int ns_cnts[num_node_sets];
    int ns_df_cnts[num_node_sets];
    int ns_ptrs[num_node_sets];
    int ns_df_ptrs[num_node_sets];
    int* ns_list       = NULL;
    double* ns_df_list = NULL;
    if (num_node_sets > 0)
    {
      new_node_sets = num_node_sets;
      if (pc_input->get_node_sets(new_node_sets, ns_ids, ns_cnts, ns_df_cnts, 
                                  ns_ptrs, ns_df_ptrs, ns_list, ns_df_list))
      {
        exo_out.put_node_sets(num_points_in, node_ids, ns_ids, 
                              ns_cnts, ns_df_cnts, ns_ptrs, ns_df_ptrs, 
                              ns_list, ns_df_list);
      }
      else
      {
        cout << "## Error: failed to write node sets!" << endl;
      }
    }
    delete [] ns_list;
    delete [] ns_df_list;
  }
  
    // read/modify/write side set if any from input file
  if (num_side_sets > 0) {
    int new_side_sets = num_side_sets;
    int ss_id_array[num_side_sets];
    int ss_cnts_array[num_side_sets];
    int ss_df_cnts_array[num_side_sets];
    int ss_ptrs[num_side_sets];
    int ss_df_ptrs[num_side_sets];
    int num_el[num_side_sets];
    int* ss_list;
    int* ss_side_list;
    int* ss_conn;
    double* ss_df_list;
    if (pc_input->get_side_sets( vol_id, new_side_sets, num_el,
                                 ss_conn, ss_id_array, ss_cnts_array,
                                 ss_df_cnts_array, ss_ptrs, ss_df_ptrs,
                                 ss_list, ss_side_list, ss_df_list ))
    {
      exo_out.put_side_sets( num_points_in, node_ids, num_el, ss_conn, 
                             num_hexes, connect, 
                             ss_id_array, ss_cnts_array, 
                             ss_df_cnts_array, ss_ptrs, ss_df_ptrs, 
                             ss_list, ss_side_list, ss_df_list );
    }
    else
    {
      cout << "## Error: failed to write side sets!" << endl;
    }
    delete [] ss_df_list;
    delete [] ss_conn;
    delete [] ss_side_list;
    delete [] ss_list;
  }

  exo_out.put_hex_blk( blk_id, num_nodes_elem, num_hexes, connect );
  exo_out.put_coor( num_points_out, x_coor, y_coor, z_coor );
}

int ReadSweepWriteSubdomains( PCExodusFile* pc_input, int vol_id, 
                              char* fileout, 
                              int num_node_sets, int num_side_sets,
                              int& num_points_out, int& num_hexes,
                              double* q_mesh, int qualityIndex, 
                              int verbose ) {
  // Read sweep subdomain parameters
  int block_id, sweep_id, num_quads, nodes_per_hex;
  pc_input->read_sweep_prop(vol_id, block_id, sweep_id, num_quads,
                            nodes_per_hex);
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
  
  int* connect = new int[num_quads * 4];
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

  // Setup CAMAL hex sweeper
  PCMLSweeper sweeper;
  sweeper.set_msg_level( verbose ? 1 : 0 );  
  sweeper.set_boundary_mesh( num_points, x_coor, y_coor, z_coor,
                             num_quads, connect,
                             num_src_surf, num_surf_quads, num_tgt_quads );
  

  // Generate swept hex mesh
  sweeper.generate_mesh( num_points_out, num_hexes );
  delete [] connect;
  connect = NULL;
  delete [] z_coor;  
  delete [] y_coor;
  delete [] x_coor;
  x_coor = y_coor = z_coor = NULL;
  
  // Retrieve mesh
  double* x_coor_m = new double[num_points_out];
  double* y_coor_m = new double[num_points_out];
  double* z_coor_m = new double[num_points_out];
  int* connect_m = new int[num_hexes * 8];
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

  // convert connectivity to exodus (PATRAN) order
  ConvertToPatranOrder( num_hexes, connect_m );

  // Write mesh
  WriteLocalExodusMesh( pc_input, vol_id, num_points, node_ids, 
                        num_points_out, num_hexes, 
                        num_node_sets, num_side_sets,
                        x_coor_m, y_coor_m, z_coor_m,
                        connect_m, sweep_id, block_id, nodes_per_hex,
                        fileout, verbose );
  delete [] connect_m;
  delete [] z_coor_m;  
  delete [] y_coor_m;
  delete [] x_coor_m;
  delete [] node_ids;

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
  time_t t0, t1;
  struct timeval tv0, tv1;

  if ( argc < 4 ) 
    {
    cout <<  "Usage: pcamal_proto <n_subdomains> <filein> <fileout>" << endl;
    return 1;
    }

  nsub = atoi( argv[1] );
  strcpy( filein, argv[2] );
  strcpy( fileout, argv[3] );

  // Start clock
  time ( &t0 );
  gettimeofday ( &tv0, NULL );

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
  int num_blks, num_node_sets, num_side_sets;
  pc_input.get_param(num_blks, num_node_sets, num_side_sets);

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
                                                num_node_sets, num_side_sets,
                                                num_points_out, num_hexes,
                                                q_mesh, PCAMAL_QUALITY_MAX_ASPECT_FROBENIUS,
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
             << " was not to be swept."
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
         <<  "  mesh quality ( "
         << qualityName.c_str()
         << " ):"
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

  // End clock
  time ( &t1 );
  gettimeofday ( &tv1, NULL );

  if ( ! myrank ) 
    {
    cout <<  endl 
	 << "# Done, walltime = " 
	 << difftime( t1, t0 )
	 << " sec ( " 
	 << difftime( tv1.tv_usec, tv0.tv_usec )
	 << " usec). Exiting PCAMAL." 
	 << endl
	 <<  "===============" 
	 << endl;
    }

  return 0;
}
