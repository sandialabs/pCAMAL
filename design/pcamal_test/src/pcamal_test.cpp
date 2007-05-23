// this is a test program for pCAMAL
// it reads a Exodus II file with pCAMAL object properties and generates
// a swept mesh for each element block in the file

#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "PCExodusFile.hpp"
#include "PCMLSweeper.hpp"

void add_high_order_nodes(int num_hexes, int* &connect, int nodes_per_hex);
void convert_hex_partran(int num_hexes, int* connect);
bool generate_filenames(int str_size, char* filein, char* fileout);
void parse_commands(int argc, char** argv, int str_size, 
                    char* filein, char* fileout, bool& verbose);

int main(int argc, char **argv)
{
  bool debug_flag = false;
  
    // parse command line arguments
  const int str_size = 128;
  char filein[str_size];
  char fileout[str_size];
  filein[0]  = 0;
  fileout[0] = 0;
  bool verbose = false;
  parse_commands(argc, argv, str_size, filein, fileout, verbose);
  
    // request input file if none entered as argument
  if (filein[0] == 0 || fileout[0] == 0)
    if (!generate_filenames(str_size, filein, fileout))
      exit (-1);
  
    // open input file
  PCExodusFile pc_input(filein, pce::read);
  int num_blks, num_node_sets, num_side_sets;
  pc_input.get_param(num_blks, num_node_sets, num_side_sets);
  if (num_blks == 0) {
    printf("Error opening/reading input file\n");
    exit(-1);
  }
  int num_vol_hexes[num_blks];
  pc_input.get_num_hexes(num_blks, num_vol_hexes);

    // for each element block 
  int vol_id;
  for (vol_id = 0; vol_id < num_blks; vol_id++) {
      // read sweep block control
    int block_id, sweep_id, num_quads;
    int nodes_per_hex;
    pc_input.read_sweep_prop(vol_id, block_id, sweep_id, num_quads,
                             nodes_per_hex);
    if (sweep_id == 0) {
        // copy the non-pCAMAL block to an output file
      continue;
    }
    printf("Sweeping pCAMAL block %d (Cubit Volume %d)\n", vol_id, sweep_id);
    
      // read coordinates (memory allocated in read_sweep_coord)
    int num_points;
    double* x_coor = NULL;
    double* y_coor = NULL;
    double* z_coor = NULL;
    int* node_ids  = NULL;
    pc_input.read_sweep_coord(vol_id, num_points, x_coor, y_coor, z_coor,
                              node_ids);
    if (num_points == 0) {
      printf("ERROR: failed to retreive coordinates for volume %d\n", vol_id);
      exit (-1);
    }
    if (debug_flag) {
      printf("\n                  ---Coordinates---\n");
      printf("  node        x            y            z\n");
      int i;
      for (i = 0; i < num_points; i++) {
        printf("%8d %12.5e %12.5e %12.5e\n", 
               i+1, x_coor[i], y_coor[i], z_coor[i]);
      }
    }

      // read connectivity
    int num_src_surf, num_lnk_surf, num_tgt_surf;
    int num_surfs = pc_input.read_sweep_surf_prop(vol_id, num_src_surf, 
                                                  num_lnk_surf, num_tgt_surf);
    if (verbose) {
      printf("\nSurface Information\n");
      printf("\tnumber sources = %d\n", num_src_surf);
      printf("\tnumber linking = %d\n", num_lnk_surf);
      printf("\tnumber target  = %d\n", num_tgt_surf);
      printf("\t         total = %d\n", num_surfs);
    }
    int *num_surf_quads = new int[num_surfs];
    pc_input.read_sweep_surf_size(vol_id, num_surfs, num_surf_quads);
    int num_tgt_quads = 0;
    int i;
    for (i = 0; i < num_src_surf; i++)
      num_tgt_quads += num_surf_quads[i];
    if (verbose) {
      printf("\nNumber of quads/surface\n");
      printf(" Surface    Quads\n");
      int i;
      for (i = 0; i < num_surfs; i++) {
        printf("%8d %8d\n", i+1, num_surf_quads[i]);
      }
    }

    int *connect = new int[num_quads * 4];
    pc_input.read_sweep_conn(vol_id, num_points, num_quads, node_ids, connect);
    if (debug_flag) {
      printf("\n           ---Connectivity---\n");
      printf("  Quad     n1        n2        n3        n4\n");
      int *c = connect;
      int i;
      for (i = 0; i < num_quads; i++) {
        printf("%8d: %8d %8d %8d %8d\n", i+1, c[0], c[1], c[2], c[3]);
        c += 4;
      }
    }  

      // setup sweeper
    int msg_level = verbose ? 1 : 0;
    PCMLSweeper sweeper;
    sweeper.set_msg_level(msg_level);
    sweeper.set_boundary_mesh(num_points, x_coor, y_coor, z_coor,
                              num_quads, connect,
                              num_src_surf, num_surf_quads, num_tgt_quads);
    delete [] num_surf_quads;
    delete [] connect;
    delete [] z_coor;
    delete [] y_coor;
    delete [] x_coor;

      // generate swept hex mesh
    int num_points_out, num_hexes;
    sweeper.generate_mesh(num_points_out, num_hexes);
    if (verbose)
      printf("Number of estimated hexes: %d -- actual: %d\n",
             num_vol_hexes[vol_id], num_hexes);

      // retrieve mesh
    x_coor = new double[num_points_out];
    y_coor = new double[num_points_out];
    z_coor = new double[num_points_out];
    int* hexes = new int[num_hexes * 8];
    sweeper.get_mesh(num_points_out, x_coor, y_coor, z_coor,
                     num_hexes, hexes);

      // convert connectivity to exodus (PATRAN) order
    convert_hex_partran(num_hexes, hexes);

      // write Exodus II output file for mesh
    int nlen = strlen(fileout);
    char *filename = new char[nlen + 10];
    sprintf(filename, "%s.vol%03d.g", fileout, sweep_id);
    PCExodusFile exo_out(filename, pce::create);
    exo_out.put_param(num_points_out, num_hexes, num_node_sets, num_side_sets);

      // read/write modified node sets
    if (num_node_sets > 0) {
      int new_node_sets = num_node_sets;
      int ns_id_array[num_node_sets];
      int ns_cnts_array[num_node_sets];
      int ns_df_cnts_array[num_node_sets];
      int ns_ptrs[num_node_sets];
      int ns_df_ptrs[num_node_sets];
      int* ns_list       = NULL;
      double* ns_df_list = NULL;
      if (!pc_input.get_node_sets(num_points, node_ids, new_node_sets, 
                                  ns_id_array, ns_cnts_array,
                                  ns_df_cnts_array, ns_ptrs, ns_df_ptrs,
                                  ns_list, ns_df_list)) {
        printf("ERROR: failed to write node sets!\n");
      }
    
      exo_out.put_node_sets(num_points, node_ids, ns_id_array, ns_cnts_array, 
                            ns_df_cnts_array, ns_ptrs, ns_df_ptrs, 
                            ns_list, ns_df_list);
      
      delete [] ns_list;
      delete [] ns_df_list;
    }
    
      // read/write modified side sets
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
      if (!pc_input.get_side_sets(vol_id, new_side_sets, num_el,
                                  ss_conn, ss_id_array, ss_cnts_array,
                                  ss_df_cnts_array, ss_ptrs, ss_df_ptrs,
                                  ss_list, ss_side_list, ss_df_list))
        printf("ERROR: failed to write side sets!\n");

      exo_out.put_side_sets(num_points, node_ids, num_el, ss_conn, 
                            num_hexes, hexes, 
                            ss_id_array, ss_cnts_array, 
                            ss_df_cnts_array, ss_ptrs, ss_df_ptrs, 
                            ss_list, ss_side_list, ss_df_list);

      delete [] ss_df_list;
      delete [] ss_conn;
      delete [] ss_side_list;
      delete [] ss_list;
    }
    delete [] node_ids;
    
      // output connectivity
    exo_out.put_hex_blk(block_id, nodes_per_hex, num_hexes, hexes);

      // output coordinates to exodus
    exo_out.put_coor(num_points_out, x_coor, y_coor, z_coor);

      // delete memory
    delete [] filename;
    delete [] hexes;
    delete [] z_coor;
    delete [] y_coor;
    delete [] x_coor;

    if (verbose)
      printf("\n");
  } // endfor each element block
}

void parse_commands(int argc, char** argv, int str_size, 
                    char* filein, char* fileout, bool& verbose)
{
  char cmd;
  int nlen = 0;
  while ((cmd = getopt(argc, argv, "hi:o:v")) != EOF) {
    switch (cmd) {
        // useage
      case 'h':
      default:
          printf("Useage: pcamal_test [OPTIONS]\n");
          printf("\t-i <filename>  input file name\n");
          printf("\t-o <filename>  output file basename\n");
          printf("\t-v             verbose output\n");
          printf("\t-h             this help message\n");          
          exit(-1);
          
            // input file name
      case 'i':
          nlen = strlen(optarg);
          strncpy(filein, optarg, str_size);
          if (nlen >= str_size)
            filein[str_size-1] = 0;
          break;

            // output file basename
      case 'o':
          nlen = strlen(optarg);
          strncpy(fileout, optarg, str_size);
          if (nlen >= str_size)
            fileout[str_size-1] = 0;
          break;

            // verbose output flag
      case 'v':
          verbose = true;
          break;
    }
  }
}

bool generate_filenames(int str_size, char* filein, char* fileout)
{
  int nlen;
  if (filein[0] == 0) {
    char filename[256];
    printf("Input filename: ");
    scanf("%s", filename);
    nlen = strlen(filename);
    if (nlen == 0) {
      printf("Error: No input file entered\n");
      return false;
    }
    strncpy(filein, filename, str_size);
    if (nlen >= str_size)
      filein[str_size-1] = 0;
  }

    // generate output file basename if none entered as argument
  if (fileout[0] == 0) {
    strncpy(fileout, filein, str_size);
  }
  return true;
}

void convert_hex_partran(int num_hexes, int* hexes)
{
  int tmp;
  int *c = hexes;
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

void add_high_order_nodes(int num_hexes, int* &hexes, int nodes_per_hex)
{
  int* new_conn = new int[num_hexes * nodes_per_hex];
  int c = 0;
  int n = 0;
  int i, j;
  for (i = 0; i < num_hexes; i++) {
    for (j = 0; j < 8; j++) {
      new_conn[n++] = hexes[c++];
    }

      // this is a test (need real node numbers)
    for (j = 8; j < nodes_per_hex; j++) {
      new_conn[n++] = -1;
    }
  }
  delete [] hexes;
  hexes = new_conn;
}

