// this is a test program for pCAMAL
// it reads a Exodus II file with pCAMAL object properties and generates
// a swept mesh for each element block in the file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "PCExodusFile.hpp"
#include "CMLSweeper.hpp"

void convert_conn(int num_hexes, int* connect);

int main(int argc, char **argv)
{
    // parse command line arguments
  char cmd;
  char *filein  = NULL;
  char *fileout = NULL;
  bool verbose = false;
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
          filein = new char[nlen+1];
          strcpy(filein, optarg);
          break;

            // output file basename
      case 'o':
          nlen = strlen(optarg);
          fileout = new char[nlen+1];
          strcpy(fileout, optarg);
          break;

            // verbose output flag
      case 'v':
          verbose = true;
          break;
    }
  }

    // request input file if none entered as argument
  if (filein == NULL) {
    char filename[256];
    printf("Input filename: ");
    scanf("%s", filename);
    nlen = strlen(filename);
    if (nlen == 0) {
      printf("Error: No input file entered\n");
      exit(-1);
    }
    filein = new char[nlen+1];
    strcpy(filein, filename);
  }

    // generate output file basename if none entered as argument
  if (fileout == NULL) {
    nlen = strlen(filein);
    fileout = new char[nlen];
    strcpy(fileout, filein);
  }

    // open input file
  PCExodusFile pc_input(filein, pce::read);
  int num_blks = pc_input.get_num_sweep_vols();

    // for each element block 
  int vol_id;
  for (vol_id = 0; vol_id < num_blks; vol_id++) {
      // read sweep block control
    int sweep_id, num_points, num_quads;
    pc_input.read_sweep_prop(vol_id, sweep_id, num_points, num_quads);
    if (sweep_id == 0)
      continue;
    if (verbose) {
      printf("\nSweeping Volume %d (Cubit Volume %d)\n", vol_id, sweep_id);
    }
    
      // read coordinates
    double *x_coor = new double[num_points];
    double *y_coor = new double[num_points];
    double *z_coor = new double[num_points];
    int *node_ids = new int[num_points];
    pc_input.read_sweep_coord(vol_id, num_points, x_coor, y_coor, z_coor,
                              node_ids);
    if (verbose) {
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
    pc_input.read_sweep_conn(vol_id, num_quads, node_ids, connect);
    if (verbose) {
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
    CMLSweeper sweeper;
    sweeper.set_boundary_mesh(num_points, x_coor, y_coor, z_coor,
                              num_quads, connect,
                              num_src_surf, num_surf_quads, num_tgt_quads);
    delete [] num_surf_quads;
    delete [] connect;
    delete [] node_ids;
    delete [] z_coor;
    delete [] y_coor;
    delete [] x_coor;

      // generate swept hex mesh
    int num_points_out, num_hexes;
    sweeper.generate_mesh(num_points_out, num_hexes);

      // retrieve mesh
    x_coor = new double[num_points_out];
    y_coor = new double[num_points_out];
    z_coor = new double[num_points_out];
    connect = new int[num_hexes * 8];
    sweeper.get_mesh(num_points_out, x_coor, y_coor, z_coor,
                     num_hexes, connect);

      // write Exodus II output file for mesh
    nlen = strlen(fileout);
    char *filename = new char[nlen + 10];
    sprintf(filename, "%s.vol%03d.g", fileout, sweep_id);
    PCExodusFile exo_out(filename, pce::create);
    exo_out.put_param(num_points_out, num_hexes);
    exo_out.put_coor(num_points_out, x_coor, y_coor, z_coor);

      // convert connectivity to exodus (PATRAN) order
    convert_conn(num_hexes, connect);
    exo_out.put_hex_blk(num_hexes, connect);
    
      // delete memory
    delete [] filename;
    delete [] connect;
    delete [] z_coor;
    delete [] y_coor;
    delete [] x_coor;
  } // endfor each element block

    // clean-up
  delete [] filein;
  delete [] fileout;
}

void convert_conn(int num_hexes, int* connect)
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
