// this is a test program for pCAMAL
// it reads a Exodus II file with pCAMAL object properties and generates
// a swept mesh for each element block in the file

#include <map>
#include <set>
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
                    char* filein, char* fileout, int& num_proc, bool& verbose);

int sweep_and_output(char* filename, int vol_id, int block_id, 
                     int nodes_per_hex, int num_node_sets, 
                     int num_side_sets, int* node_ids, 
                     PCExodusFile& input, PCMLSweeper& sweeper);

bool sub_vol_input(int num_points, double* x_coor, double* y_coor,
                   double* z_coor, int num_quads_out, int* connect,
                   std::vector<int>& src_ids, std::vector<int>& lnk_ids, 
                   std::vector<int>& tgt_ids, int& num_points_sub, 
                   double* &x_sub, double* &y_sub, double* &z_sub,
                   int& num_src_surfs_sub, int* &connect_sub,
                   int* &num_surf_quads_sub, int& num_tgt_quads_sub);

void xform_node_set(int num_node_sets, int num_nodes, int* node_ids,
                    PCExodusFile& in, PCExodusFile& out);

void xform_side_set(int index, int num_side_sets, int num_points,
                    int* node_ids, int num_hexes, int* hexes,
                    PCExodusFile& in, PCExodusFile& out);

int main(int argc, char **argv)
{
  bool debug_flag = false;
  
    // parse command line arguments
  const int str_size = 128;
  char filein[str_size];
  char fileout[str_size];
  filein[0]  = 0;
  fileout[0] = 0;
  int  num_proc = 1;
  bool verbose = false;
  parse_commands(argc, argv, str_size, filein, fileout, num_proc, verbose);
  
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

    // compute the average number of hexes per processor
  int num_vol_hexes[num_blks];
  int total_hexes = pc_input.get_num_hexes(num_blks, num_vol_hexes);
  int ave_hexes = 
      num_proc > 1 ? (int)(total_hexes/num_proc + 0.5) : total_hexes;

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
//     int msg_level = verbose ? 1 : 0;
    int msg_level = 0;
    PCMLSweeper sweeper(NULL, msg_level);
    sweeper.set_boundary_mesh(num_points, x_coor, y_coor, z_coor,
                              num_quads, connect,
                              num_src_surf, num_surf_quads, num_tgt_quads);
    delete [] num_surf_quads;
    delete [] connect;
    delete [] z_coor;
    delete [] y_coor;
    delete [] x_coor;
    x_coor = y_coor = z_coor = NULL;
    num_surf_quads = connect = NULL;

      // load balancing here
    int num_sub_vol = (int)(num_vol_hexes[vol_id]/ave_hexes + 0.5);
    if (num_sub_vol > 1) {
      int num_points_out = 0;
      int num_quads_out = 0;
      bool ret_val = sweeper.balance_mesh(num_sub_vol, num_points_out, 
                                          num_quads_out);

        // retrieve load balanced coordinates and connectivity
      if (ret_val) {
        x_coor = new double[num_points_out];
        y_coor = new double[num_points_out];
        z_coor = new double[num_points_out];
        connect = new int[num_quads_out * 4];
        ret_val = sweeper.get_shell(num_points_out, x_coor, y_coor,
                                    z_coor, num_quads_out, connect);
      }

        // generate a hex mesh for each sub-volume
      if (ret_val) {
        int j;
        for (j = 0; j < num_sub_vol; j++) {
            // retrieve the source, linking, and target quad ids
          std::vector<int> src_ids;
          std::vector<int> lnk_ids;
          std::vector<int> tgt_ids;
          int num_layers;
          ret_val = sweeper.get_block_quads(j, num_layers, src_ids, lnk_ids,
                                            tgt_ids);

            // prepare individual sub-volume data and mesh it
          int num_points_sub;
          int num_quads_sub = src_ids.size() + lnk_ids.size() + tgt_ids.size();
          int num_src_surfs_sub;
          int num_tgt_quads_sub;
          int* connect_sub = NULL;
          int* num_surf_quads_sub = NULL;
          double* x_sub = NULL;
          double* y_sub = NULL;
          double* z_sub = NULL;
          if (ret_val) {
            ret_val = sub_vol_input(num_points_out, x_coor, y_coor,
                                    z_coor, num_quads_out, connect,
                                    src_ids, lnk_ids, tgt_ids,
                                    num_points_sub, x_sub, y_sub, z_sub,
                                    num_src_surfs_sub, connect_sub,
                                    num_surf_quads_sub, num_tgt_quads_sub);
          }
          
          if (ret_val) {
            PCMLSweeper sweep2(NULL, msg_level);
            sweep2.set_boundary_mesh(num_points_sub, x_sub, y_sub, z_sub,
                                     num_quads_sub, connect_sub,
                                     num_src_surfs_sub, num_surf_quads_sub,
                                     num_tgt_quads_sub);

            char filename[256];
            sprintf(filename, "%s.vol%03d-%02d.g", fileout, sweep_id, j+1);

            int num_hexes = sweep_and_output(filename, vol_id, block_id, 
                                             nodes_per_hex, num_node_sets, 
                                             num_side_sets, node_ids, 
                                             pc_input, sweep2);
            if (verbose && num_hexes > 0)
              printf("Number of estimated hexes: %d: sub%02d %d\n",
                     num_vol_hexes[vol_id], j+1, num_hexes);
          }

            // delete memory
          delete [] z_sub;
          delete [] y_sub;
          delete [] x_sub;
          delete [] num_surf_quads_sub;
          delete [] connect_sub;
          x_sub = y_sub = z_sub = NULL;
          num_surf_quads_sub = connect_sub = NULL;
        }
      }     
        // clean-up more memory
      delete [] connect;
      delete [] z_coor;
      delete [] y_coor;
      delete [] x_coor;
      x_coor = y_coor = z_coor = NULL;
      connect = NULL;
    }
    
    else {
      char filename[256];
      sprintf(filename, "%s.vol%03d.g", fileout, sweep_id);

      int num_hexes = sweep_and_output(filename, vol_id, block_id, 
                                       nodes_per_hex, num_node_sets, 
                                       num_side_sets, node_ids, 
                                       pc_input, sweeper);
      if (verbose && num_hexes > 0)
        printf("Number of estimated hexes: %d -- actual: %d\n",
               num_vol_hexes[vol_id], num_hexes);
    }
    if (verbose)
      printf("\n");

      // clean-up memory
    delete [] node_ids;    
  } // endfor each element block
}

void parse_commands(int argc, char** argv, int str_size, 
                    char* filein, char* fileout, int& num_proc,
                    bool& verbose)
{
  char cmd;
  int nlen = 0;
  while ((cmd = getopt(argc, argv, "hi:n:o:v")) != EOF) {
    switch (cmd) {
        // useage
      case 'h':
      default:
          printf("Useage: pcamal_test [OPTIONS]\n");
          printf("\t-i <filename>  input file name\n");
          printf("\t-o <filename>  output file basename\n");
          printf("\t-n <number>    use number of processors\n");
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

            // number of processors
      case 'n':
          num_proc = atoi(optarg);
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

void xform_node_set(int num_node_sets, int num_points, int* node_ids,
                    PCExodusFile& input, PCExodusFile& out)
{
  int new_node_sets = num_node_sets;
  int ns_id_array[num_node_sets];
  int ns_cnts_array[num_node_sets];
  int ns_df_cnts_array[num_node_sets];
  int ns_ptrs[num_node_sets];
  int ns_df_ptrs[num_node_sets];
  int* ns_list       = NULL;
  double* ns_df_list = NULL;
  if (input.get_node_sets(new_node_sets, ns_id_array,
                          ns_cnts_array, ns_df_cnts_array,
                          ns_ptrs, ns_df_ptrs, ns_list, ns_df_list)) {
    out.put_node_sets(num_points, node_ids, ns_id_array, 
                      ns_cnts_array, ns_df_cnts_array, 
                      ns_ptrs, ns_df_ptrs, ns_list, ns_df_list);
  }
  else {
    printf("ERROR: failed to write node sets!\n");
  }
  delete [] ns_list;
  delete [] ns_df_list;
}

void xform_side_set(int index, int num_side_sets, int num_points, 
                    int* node_ids, int num_hexes, int* hexes,
                    PCExodusFile& input, PCExodusFile& out)
{
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
  if (input.get_side_sets(index, new_side_sets, num_el,
                          ss_conn, ss_id_array, ss_cnts_array,
                          ss_df_cnts_array, ss_ptrs, ss_df_ptrs,
                          ss_list, ss_side_list, ss_df_list)) {
    out.put_side_sets(num_points, node_ids, num_el, ss_conn, 
                      num_hexes, hexes, 
                      ss_id_array, ss_cnts_array, 
                      ss_df_cnts_array, ss_ptrs, ss_df_ptrs, 
                      ss_list, ss_side_list, ss_df_list);
  }
  else {
    printf("ERROR: failed to write side sets!\n");
  }
  delete [] ss_df_list;
  delete [] ss_conn;
  delete [] ss_side_list;
  delete [] ss_list;
}

int sweep_and_output(char* filename, int vol_id, int block_id, 
                     int nodes_per_hex, int num_node_sets, int num_side_sets,
                     int* node_ids, PCExodusFile& input, PCMLSweeper& sweeper)
{
    // generate swept hex mesh
  int num_points_out;
  int num_hexes = 0;
  bool ret_val = sweeper.generate_mesh(num_points_out, num_hexes);

  double x_coor[num_points_out];
  double y_coor[num_points_out];
  double z_coor[num_points_out];
  int    hexes[num_hexes * 8];

    // retrieve mesh
  if (ret_val) {
    ret_val = sweeper.get_mesh(num_points_out, x_coor, y_coor, z_coor,
                               num_hexes, hexes);
  }

    // convert connectivity to exodus (PATRAN) order
  if (ret_val) {
    convert_hex_partran(num_hexes, hexes);

      // write Exodus II output file for mesh
    PCExodusFile exo_out(filename, pce::create);
    int err = exo_out.put_param(num_points_out, num_hexes, num_node_sets, 
                                num_side_sets);

      // read/write modified node sets
    if (err == 0 && num_node_sets > 0) {
      xform_node_set(num_node_sets, num_points_out, node_ids, input, exo_out);
    }
    
      // read/write modified side sets
    if (err == 0 && num_side_sets > 0) {
      xform_side_set(vol_id, num_side_sets, num_points_out, node_ids,
                     num_hexes, hexes, input, exo_out);
    }
    
      // output connectivity
    if (err == 0)
      err = exo_out.put_hex_blk(block_id, nodes_per_hex, num_hexes, hexes);
  
      // output coordinates to exodus
    if (err == 0)
      err = exo_out.put_coor(num_points_out, x_coor, y_coor, z_coor);

    num_hexes = err == 0 ? num_hexes : 0;
  }

  return num_hexes;
}

bool sub_vol_input(int num_points, double* x_coor, double* y_coor,
                   double* z_coor, int num_quads, int* connect,
                   std::vector<int>& src_ids, std::vector<int>& lnk_ids, 
                   std::vector<int>& tgt_ids, int& num_points_out, 
                   double* &x, double* &y, double* &z,
                   int& num_src_surfs, int* &quads,
                   int* &num_surf_quads, int& num_tgt_quads)
{
    // compute total number of quads in sub-volume
  assert(src_ids.size() == tgt_ids.size());
  int num_quads_out = src_ids.size() + lnk_ids.size() + tgt_ids.size();
  if (num_quads_out <= 0)
    return false;
  
    // compute number of quad on souce, linking, target surfaces
  num_src_surfs = 1;
  num_surf_quads = new int[3];
  num_surf_quads[0] = src_ids.size();
  num_surf_quads[1] = lnk_ids.size();
  num_surf_quads[2] = tgt_ids.size();
  num_tgt_quads = num_surf_quads[2];
  
    // copy quads first
  std::set<int> used_ids;
  quads = new int[num_quads_out * 4];
  int cnt = 0;
  int i, j;
  for (i = 0; i < src_ids.size(); i++) {
    int index = (src_ids[i] - 1) * 4;
    for (j = 0; j < 4; j++) {
      int id = connect[index + j];
      quads[cnt++] = id;
      used_ids.insert(id);
    }
  }
  for (i = 0; i < lnk_ids.size(); i++) {
    int index = (lnk_ids[i] - 1) * 4;
    for (j = 0; j < 4; j++) {
      int id = connect[index + j];
      quads[cnt++] = id;
      used_ids.insert(id);
    }
  }
  for (i = 0; i < tgt_ids.size(); i++) {
    int index = (tgt_ids[i] - 1) * 4;
    for (j = 0; j < 4; j++) {
      int id = connect[index + j];
      quads[cnt++] = id;
      used_ids.insert(id);
    }
  }
  assert(cnt/4 == num_quads_out);
  
  num_points_out = used_ids.size();
  if (num_points_out <= 0)
    return false;

    // copy the coordinates
  std::map<int,int> node_map;
  cnt = 0;
  x = new double[num_points_out];
  y = new double[num_points_out];
  z = new double[num_points_out];
  std::set<int>::iterator it;
  for (it = used_ids.begin(); it != used_ids.end(); it++) {
    int index = (*it);
    node_map[index] = cnt;
    x[cnt] = x_coor[index];
    y[cnt] = y_coor[index];
    z[cnt] = z_coor[index];
    cnt++;
  }
  assert(cnt == num_points_out);

    // convert connectivity to new node numbering
  for (i = 0; i < num_quads_out * 4; i++) {
    quads[i] = node_map[quads[i]];
  }
  
  return true;
}
