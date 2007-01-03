// this program reads a file generated by Cubit with element blocks,
// nodesets and sidesets for sweepable volumes (defined by the user)
// and generates a pcamal format Exodus II file.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "exodusII.h"

// forward declarations
int copy_init_rec(int exo_in, int exo_out, bool update,
                  int &num_nodes, int &num_elem, int &num_blks,
                  int &num_node_sets, int &num_side_sets);
int copy_info_rec(int exo_in, int exo_out, bool update);
int copy_qa_rec(int exo_in, int exo_out, bool update);
int copy_coord_rec(int exo_in, int exo_out, int num_nodes, 
                   bool update, bool verbose);
int copy_node_sets(int exo_in, int exo_out, int num_nodes_sets, 
                   bool update, bool verbose);
int copy_elem_blk_rec(int exo_in, int exo_out, int num_blks,
                      bool update, bool verbose);
int copy_side_sets(int exo_in, int exo_out, int num_side_sets, 
                   bool update, bool verbose);
int copy_elem_order_map(int exo_in, int exo_out, int num_quads, 
                        bool update, bool verbose);
int copy_obj_prop_rec(int exo_in, int exo_out, int num_blks, 
                      int num_node_sets, int num_side_sets, 
                      bool update, bool verbose);
void print_property(int num_ss_props, char *ss_prop_name, int*ss_prop_values);


int main(int argc, char **argv)
{
  bool verbose = false;
  bool update  = false;
  
  char cmd;
  char *filein  = NULL;
  char *fileout = NULL;
  while ( (cmd = getopt(argc, argv, "f:ho:uv")) != EOF) {
    switch (cmd) {
      case 'v':
          verbose = true;
          break;
      case 'u':
          update = true;
          break;
      case 'f':
          filein = optarg;
          break;
      case 'o':
          fileout = optarg;
          break;
      case 'h':
      default:
          printf("Usage: bisected-file [OPTION]\n");
          printf("   -f <filename>  read from this file\n");
          printf("   -o <filename>  write to this file\n");
          printf("   -u             update file QA and object properties\n");
          printf("   -v             verbose output\n");
          printf("   -h             this help message\n");
          exit(-1);
    }
  }
          

  char input_name[MAX_LINE_LENGTH+1];
  if (filein == NULL) {
      printf("Cubit input filename: ");
      scanf("%s", input_name);
    }
  else {
    strcpy(input_name, filein);
  }

  char output_name[MAX_LINE_LENGTH+1];
  if (fileout == NULL) {
    if (update) {
      printf("pCAMAL output filename: ");
      scanf("%s", output_name);
    }
  }
  else {
    strcpy(output_name, fileout);
    update = true;
  }
  
    // open file for reading
  float version;
  int cpu_word = 8;
  int file_word = 0;
  int exo_in = ex_open(input_name, EX_READ, &cpu_word, &file_word, &version);
  if (exo_in < 0) {
    printf("Error opening Exodus2 file: %s\n", input_name);
    exit(-1);
  }
  else {
    printf("Reading from file: %s\n", input_name);
  }

    // open file for writing
  int exo_out = 0;
  if (update) {
    exo_out = ex_create(output_name, EX_CLOBBER, &cpu_word, &file_word);
    if (exo_out < 0) {
      printf("Error opening Exodus2 file: %s\n", output_name);
      exit(-1);
    }
    else {
      printf("Writing to file: %s\n", output_name);
    }
  }

    // read/write header
  int num_nodes, num_elems, num_blks, num_node_sets, num_side_sets;
  int error = copy_init_rec(exo_in, exo_out, update, num_nodes, num_elems,
                            num_blks, num_node_sets, num_side_sets);
  
    // read/write info records
  if (error == 0)
    error = copy_info_rec(exo_in, exo_out, update);
  
    // read/write QA records
  if (error == 0)
    error = copy_qa_rec(exo_in, exo_out, update);
  
    // read/write coordinate names
  if (error == 0)
    error = copy_coord_rec(exo_in, exo_out, num_nodes, update, verbose);

    // read/write node sets
  if (error == 0)
    error = copy_node_sets(exo_in, exo_out, num_node_sets, update, verbose);

    // read/write element blocks
  if (error == 0)
    error = copy_elem_blk_rec(exo_in, exo_out, num_blks, update, verbose);

    // read/write side sets
  if (error == 0)
    error = copy_side_sets(exo_in, exo_out, num_side_sets, update, verbose);

    // read/write element order map
  if (error == 0)
    error = copy_elem_order_map(exo_in, exo_out, num_elems, update, verbose);
  
    // write object properties for sweepable blocks
  if (error == 0)
    error = copy_obj_prop_rec(exo_in, exo_out, num_blks, num_node_sets,
                              num_side_sets, update, verbose);
  
    // close file
  if (ex_close(exo_in) != 0)
    printf("Error closing input file\n");
  if (update && ex_close(exo_out) != 0)
    printf("Error closing output file\n");

  exit(error);
}

int copy_init_rec(int exo_in, int exo_out, bool update,
                  int &num_nodes, int &num_elem, int &num_blks,
                  int &num_node_sets, int &num_side_sets)
{
  char title[MAX_LINE_LENGTH+1];
  int num_dim;
  int error = ex_get_init(exo_in, title, &num_dim, &num_nodes, &num_elem,
                          &num_blks, &num_node_sets, &num_side_sets);
  if (error == 0) {
    printf("\nTitle: %s\n", title);
    printf("\tDimensions:           %d\n", num_dim);
    printf("\tNum. of nodes:        %d\n", num_nodes);
    printf("\tNum. of elements:     %d\n", num_elem);
    printf("\tNum. of elem. blocks: %d\n", num_blks);
    printf("\tNum. of node sets:    %d\n", num_node_sets);
    printf("\tNum. of side sets:    %d\n", num_side_sets);    
  }
  else {
    printf("\nError reading initialization parameters\n");
    exit(error);
  }

    // write header
  if (update > 0 && error == 0) {
    int len = strlen(title);
    if (len + 20 < MAX_LINE_LENGTH + 1)
      strcat(title, " -- pCAMAL interface");
    error = ex_put_init(exo_out, title, num_dim, num_nodes, num_elem,
                        num_blks, num_node_sets, num_side_sets);
  }
  return error;
}

int copy_info_rec(int exo_in, int exo_out, bool update)
{
  int num_info;
  float fdum;
  char* cdum;
  int error = ex_inquire(exo_in, EX_INQ_INFO, &num_info, &fdum, cdum);
  if (error == 0) {
    if (num_info > 0) {
      char **info = new char*[num_info];
      int i;
      for (i = 0; i < num_info; i++)
        info[i] = new char[MAX_LINE_LENGTH+1];
      error = ex_get_info(exo_in, info);
      if (error == 0) {
        printf("\n");
        for (i = 0; i < num_info; i++)
          printf("Info record %d: %s\n", i+1, info[i]);
      }
      else {
        printf("\nError reading info records\n");
      }

        // copy info to output
      if (update && error == 0)
        error = ex_put_info(exo_out, num_info, info);
      
      for (i = 0; i < num_info; i++)
        delete [] info[i];
      delete [] info;
    }
    else {
      printf("\nNo info records in this file.\n");
    }
  }
  else {
    printf("\nError inquiring about info records\n");
  }
  return error;
}

int copy_qa_rec(int exo_in, int exo_out, bool update)
{
  int num_qa;
  float fdum;
  char *cdum;
  int error = ex_inquire(exo_in, EX_INQ_QA, &num_qa, &fdum, cdum);
  if (error == 0) {
    if (num_qa > 0) {
      const int MAX_QA_REC = 10;
      if (num_qa > MAX_QA_REC) {
        printf("\nToo many QA records (%d); will list only %d\n", 
               num_qa, MAX_QA_REC);
        num_qa = MAX_QA_REC;
      }
      else {
        printf("\n");
      }
        
      char *qa[MAX_QA_REC][4];
      int i, j;
      for (i = 0; i <= num_qa; i++)
        for (j = 0; j < 4; j++)
          qa[i][j] = new char[MAX_STR_LENGTH+1];
      error = ex_get_qa(exo_in, qa);
      if (error == 0) {
        for (i = 0; i < num_qa; i++) {
          printf("QA record %d\n", i+1);
          printf("\tAnalysis code:      %s\n", qa[i][0]);
          printf("\tCode QA Descriptor: %s\n", qa[i][1]);
          printf("\tAnalysis date:      %s\n", qa[i][2]);
          printf("\tAnalysis time:      %s\n", qa[i][3]);
        }

        if (update && num_qa+1 < MAX_QA_REC) {
            // get date and time
          time_t result = time(NULL);
          struct tm *timeptr = localtime(&result);
          char camal_date[MAX_STR_LENGTH+1];
          strftime(camal_date, MAX_STR_LENGTH, "%D", timeptr);
          char camal_time[MAX_STR_LENGTH+1];
          strftime(camal_time, MAX_STR_LENGTH, "%H:%M:%S", timeptr);
    
            // write pCamal QA record
          strcpy(qa[num_qa][0], "BISECTED-FILE");
          strcpy(qa[num_qa][1], "1.2");
          strcpy(qa[num_qa][2], camal_date);
          strcpy(qa[num_qa][3], camal_time);
          ++num_qa;
          error = ex_put_qa(exo_out, num_qa, qa);
          if (error != 0)
            printf("Error writing pCAMAL QA record\n");
        }
      }
      else {
        printf("Error reading QA records\n");
      }
      for (i = 0; i < num_qa; i++)
        for (j = 0; j < 4; j++)
          delete [] qa[i][j];
    }
    else {
      printf("No QA records in this file.\n");
    }
  }
  else {
    printf("\nError inquiring about QA records\n");
  }
  return error;
}

int copy_coord_rec(int exo_in, int exo_out, int num_nodes, 
                   bool update, bool verbose)
{
  char *coord_names[3];
  int i;
  for (i = 0; i < 3; i++)
    coord_names[i] = new char[MAX_STR_LENGTH+1];
  int error = ex_get_coord_names(exo_in, coord_names);
  if (error == 0) {
    printf("\n node %12s %12s %12s\n", 
           coord_names[0], coord_names[1], coord_names[2]);
  }
  else {
    printf("\n node      x-coord      y-coord      z-coord\n");
  }

    // write coordinate names
  if (error == 0 && update)
    error = ex_put_coord_names(exo_out, coord_names);
  
  for (i = 0; i < 3; i++)
    delete [] coord_names[i];
  
    // read coordinates
  double *x_coor = new double[num_nodes];
  double *y_coor = new double[num_nodes];
  double *z_coor = new double[num_nodes];
  error = ex_get_coord(exo_in, x_coor, y_coor, z_coor);
  if (error == 0) {
    int i;
    for (i = 0; i < num_nodes; i++) {
      if (!verbose && i > 99)
        break;
      printf("%5d %12.5f %12.5f %12.5f\n", 
             i+1, x_coor[i], y_coor[i], z_coor[i]);
    }
  }
  else {
    printf("Error reading coordinates\n");
  }

    // write coordinates
  if (error == 0 && update)
    error = ex_put_coord(exo_out, x_coor, y_coor, z_coor);

  delete [] x_coor;
  delete [] y_coor;
  delete [] z_coor;

  return error;
}

int copy_elem_blk_rec(int exo_in, int exo_out, int num_blks, 
                      bool update, bool verbose)
{
  int *eb_ids = new int[num_blks];
  int error = ex_get_elem_blk_ids(exo_in, eb_ids);
  if (error == 0) {
      // read the blocks
    int i;
    for (i = 0; i < num_blks; i++) {
      int num_elem_blk, num_node_elem, num_attr;
      char elem_type[MAX_STR_LENGTH+1];
      error = ex_get_elem_block(exo_in, eb_ids[i], elem_type, &num_elem_blk,
                                &num_node_elem, &num_attr);

        // write element block parameters
      if (error == 0 && update)
        error = ex_put_elem_block(exo_out, eb_ids[i], elem_type, num_elem_blk,
                                  num_node_elem, num_attr);

      if (error == 0) {
        printf("\nElement block id:    %d\n", eb_ids[i]);
        printf("\tElement type:        %s\n", elem_type);
        printf("\tNum. elem. in block: %d\n", num_elem_blk);
        printf("\tNodes per elem.:     %d\n", num_node_elem);
        printf("\tNum. attributes:     %d\n", num_attr);

        int *connect = new int[num_node_elem * num_elem_blk];
        error = ex_get_elem_conn(exo_in, eb_ids[i], connect);
        if (error == 0) {
          printf("\n    elem        n1       n2       n3       n4\n");
          int *c = connect;
          int j;
          for (j = 0; j < num_elem_blk; j++) {
            if (!verbose && j > 99)
              break;
            printf("%8d: %8d %8d %8d %8d\n", j+1, c[0], c[1], c[2], c[3]);
            c += 4;
          }
        }
        else {
          printf("Error reading element connectivity\n");
        }

          // write element block connectivities
        if (error == 0 && update)
          error = ex_put_elem_conn(exo_out, eb_ids[i], connect);
        
        delete [] connect;

          // read element block attributes
        double *attrib = NULL;
        if (error == 0) {
          attrib = new double[num_attr * num_elem_blk];
          error = ex_get_elem_attr(exo_in, eb_ids[i], attrib);
        }
          // write element block attributes
        if (error == 0 && update) 
          error = ex_put_elem_attr(exo_out, eb_ids[i], attrib);
          
        delete [] attrib;
      }
      else {
        printf("Error reading element block\n");
      }
    }
  }
  else {
    printf("Error reading element block ids\n");
  }
  delete [] eb_ids;

  return error;
}

int copy_obj_prop_rec(int exo_in, int exo_out, int num_blks, 
                      int num_node_sets, int num_side_sets,
                      bool update, bool verbose)
{
  float fdum;
  char* cdum = NULL;
  int num_eb_props = 0;
  int num_ns_props = 0;
  int num_ss_props = 0;
  int error = ex_inquire(exo_in, EX_INQ_EB_PROP, &num_eb_props, &fdum, cdum);
  if (error == 0)
    error = ex_inquire(exo_in, EX_INQ_NS_PROP, &num_ns_props, &fdum, cdum);
  if (error == 0)
    error = ex_inquire(exo_in, EX_INQ_SS_PROP, &num_ss_props, &fdum, cdum);
  
  if (error != 0) {
    printf("Error reading number of object properties\n");
    return error;
  }
  
    // allocate memory for object properties
  int i;
  char** eb_prop_names = new char*[num_eb_props];
  for (i = 0; i < num_eb_props; i++)
    eb_prop_names[i] = new char[MAX_STR_LENGTH + 1];
  char** ns_prop_names = new char*[num_ns_props];
  for (i = 0; i < num_ns_props; i++)
    ns_prop_names[i] = new char[MAX_STR_LENGTH + 1];
  char** ss_prop_names = new char*[num_ss_props];
  for (i = 0; i < num_ss_props; i++)
    ss_prop_names[i] = new char[MAX_STR_LENGTH + 1];
  int** eb_prop_values = new int*[num_eb_props];
  for (i = 0; i < num_eb_props; i++)
    eb_prop_values[i] = new int[num_blks];
  int** ns_prop_values = new int*[num_ns_props];
  for (i = 0; i < num_ns_props; i++)
    ns_prop_values[i] = new int[num_node_sets];
  int** ss_prop_values = new int*[num_ss_props];
  for (i = 0; i < num_ss_props; i++)
    ss_prop_values[i] = new int[num_side_sets];

    // read the properties
  if (error == 0)
    error = ex_get_prop_names(exo_in, EX_ELEM_BLOCK, eb_prop_names);
  for (i = 0; i < num_eb_props && error == 0; i++)
    error = ex_get_prop_array(exo_in, EX_ELEM_BLOCK, eb_prop_names[i],
                              eb_prop_values[i]);
  if (error == 0)
    error = ex_get_prop_names(exo_in, EX_NODE_SET, ns_prop_names);
  for (i = 0; i < num_ns_props && error == 0; i++)
    error = ex_get_prop_array(exo_in, EX_NODE_SET, ns_prop_names[i],
                              ns_prop_values[i]);
  if (error == 0)
    error = ex_get_prop_names(exo_in, EX_SIDE_SET, ss_prop_names);
  for (i = 0; i < num_ss_props && error == 0; i++)
    error = ex_get_prop_array(exo_in, EX_SIDE_SET, ss_prop_names[i],
                              ss_prop_values[i]);
    

    // write the object parameters
  if (update && error == 0) {
      // element block properties
    if (num_eb_props > 1) {
      if (error == 0)
        error = ex_put_prop_names(exo_out, EX_ELEM_BLOCK, num_eb_props - 1, 
                                  &eb_prop_names[1]);

      if (verbose && error == 0)
        printf("\nNo. of block properties:   %d\n", num_eb_props);

      for (i = 0; i < num_eb_props && error == 0; i++) {
        error = ex_put_prop_array(exo_out, EX_ELEM_BLOCK, eb_prop_names[i], 
                                  eb_prop_values[i]);
        if (verbose) {
          print_property(num_blks, eb_prop_names[i], eb_prop_values[i]);
        }
      }
    }
  
      // node set properties
    if (num_ns_props > 1) {
      if (error == 0)
        error = ex_put_prop_names(exo_out, EX_NODE_SET, num_ns_props - 1, 
                                  &ns_prop_names[1]);

      if (verbose && error == 0) 
        printf("\nNo. of nodeset properties: %d\n", num_ns_props);

      for (i = 0; i < num_ns_props && error == 0; i++) {
        error = ex_put_prop_array(exo_out, EX_NODE_SET, ns_prop_names[i], 
                                  ns_prop_values[i]);
        if (verbose) {
          print_property(num_node_sets, ns_prop_names[i], ns_prop_values[i]);
        }
      }
    }

      // side set properties
    if (num_ss_props > 1) {
      if (error == 0)
        error = ex_put_prop_names(exo_out, EX_SIDE_SET, num_ss_props - 1,
                                  &ss_prop_names[1]);

      if (verbose && error == 0)
        printf("\nNo. of sideset properties: %d\n", num_ss_props);

      for (i = 0; i < num_ss_props && error == 0; i++) {
        error = ex_put_prop_array(exo_out, EX_SIDE_SET, ss_prop_names[i], 
                                  ss_prop_values[i]);
        if (verbose) {
          print_property(num_side_sets, ss_prop_names[i], ss_prop_values[i]);
        }
      }
    }
  }
  
    // delete allocate memory
  for (i = 0; i < num_eb_props; i++)
    delete [] eb_prop_names[i];
  delete [] eb_prop_names;
  for (i = 0; i < num_ns_props; i++)
    delete [] ns_prop_names[i];
  delete [] ns_prop_names;
  for (i = 0; i < num_ss_props; i++)
    delete [] ss_prop_names[i];
  delete [] ss_prop_names;
  for (i = 0; i < num_eb_props; i++)
    delete [] eb_prop_values[i];
  delete [] eb_prop_values;
  for (i = 0; i < num_ns_props; i++)
    delete [] ns_prop_values[i];
  delete [] ns_prop_values;
  for (i = 0; i < num_ss_props; i++)
    delete [] ss_prop_values[i];
  delete [] ss_prop_values;

  return error;
}


int copy_node_sets(int exo_in, int exo_out, int num_node_sets, 
                   bool update, bool verbose)
{
  int error = 0;
  int *ns_ids = NULL;
  if (num_node_sets > 0) {
    ns_ids = new int[num_node_sets];
    error = ex_get_node_set_ids(exo_in, ns_ids);
    int i;
    for (i = 0; i < num_node_sets && error == 0; i++) {
        // read node set parameters
      int num_nodes_in_set, num_dist_in_set;
      error = ex_get_node_set_param(exo_in, ns_ids[i], &num_nodes_in_set,
                                    &num_dist_in_set);
      if (error == 0) {
        printf("\nNodes set %d:\n", ns_ids[i]);
        printf("\tnumber nodes: %d\n", num_nodes_in_set);
        printf("\tnumber df's:  %d\n", num_dist_in_set);
      }

        // write nodes_set parameters
      if (error == 0 && update)
        error = ex_put_node_set_param(exo_out, ns_ids[i], num_nodes_in_set,
                                      num_dist_in_set);
      
        // read node set
      int *node_set_list = new int[num_nodes_in_set];
      if (error == 0) {
        error = ex_get_node_set(exo_in, ns_ids[i], node_set_list); 
      }
      if (error == 0) {
        printf("\n");
        int cnt = 0;
        int j;
        for (j = 0; j < num_nodes_in_set; j++) {
          if (!verbose && j > 99)
            break;
          
          if (cnt == 0) {
            printf("\n%5d:    ", j + 1);
            cnt = 10;
          }
          printf(" %5d", node_set_list[j]);
          --cnt;
        }
        printf("\n");
      }

        // write node set
      if (error == 0 && update) {
        error = ex_put_node_set(exo_out, ns_ids[i], node_set_list);
      }
      delete [] node_set_list;

        // distribution factors
      if (num_dist_in_set > 0) {
        
          // read node set distribution factors
        double *df_list = NULL;
        if (error == 0) {
          df_list = new double[num_nodes_in_set];
          error = ex_get_node_set_dist_fact(exo_in, ns_ids[i], df_list);
        }

          // write node set distribution factors
        if (error == 0 && update) 
          error = ex_put_node_set_dist_fact(exo_out, ns_ids[i], df_list);

        delete [] df_list;
      }
    }
  }
  
    // clean up
  delete [] ns_ids;

  return error;
}

int copy_side_sets(int exo_in, int exo_out, int num_side_sets, 
                   bool update, bool verbose)
{
  int error = 0;
  int *ss_ids = NULL;
  if (num_side_sets > 0) {
    ss_ids = new int[num_side_sets];
    error = ex_get_side_set_ids(exo_in, ss_ids);
    int i;
    for (i = 0; i < num_side_sets && error == 0; i++) {
        // read side set parameters
      int num_sides_in_set, num_dist_in_set;
      error = ex_get_side_set_param(exo_in, ss_ids[i], &num_sides_in_set,
                                    &num_dist_in_set);
      if (error == 0) {
        printf("\nSides set %d:\n", ss_ids[i]);
        printf("\tnumber sides: %d\n", num_sides_in_set);
        printf("\tnumber df's:  %d\n", num_dist_in_set);
      }

        // write sides_set parameters
      if (error == 0 && update)
        error = ex_put_side_set_param(exo_out, ss_ids[i], num_sides_in_set,
                                      num_dist_in_set);
      
        // read side set
      int *elem_list = new int[num_sides_in_set];
      int *side_list = new int[num_sides_in_set];
      if (error == 0) {
        error = ex_get_side_set(exo_in, ss_ids[i], elem_list, side_list); 
      }

        // second block need special processing for bisected_cyl
      const int OFFSET = 1942;
      const int BISECTED_ID = 110;
      if (error == 0 && ss_ids[i] == BISECTED_ID) {
        int j;
        for (j = 0; j < num_sides_in_set; j++) {
          if (elem_list[j] > OFFSET)
            elem_list[j] -= OFFSET;
        }
      }
      
      if (error == 0) {
        printf("\n");
        int cnt = 0;
        int j;
        for (j = 0; j < num_sides_in_set; j++) {
          if (!verbose && j > 99)
            break;
          
          if (cnt == 0) {
            printf("\n%5d:    ", j + 1);
            cnt = 5;
          }
          printf(" (%4d:%4d)", elem_list[j], side_list[j]);
          --cnt;
        }
        printf("\n");
      }

        // write side set
      if (error == 0 && update) {
        error = ex_put_side_set(exo_out, ss_ids[i], elem_list, side_list);
      }
      delete [] side_list;
      delete [] elem_list;

        // distribution factors
      if (num_dist_in_set > 0) {
        
          // read side set distribution factors
        double *df_list = NULL;
        if (error == 0) {
          df_list = new double[num_dist_in_set];
          error = ex_get_side_set_dist_fact(exo_in, ss_ids[i], df_list);
        }

          // write side set distribution factors
        if (error == 0 && update) 
          error = ex_put_side_set_dist_fact(exo_out, ss_ids[i], df_list);

        delete [] df_list;
      }
    }
  }
  
    // clean up
  delete [] ss_ids;
  
  return error;
}

int copy_elem_order_map(int exo_in, int exo_out, int num_quads, 
                        bool update, bool verbose)
{
  int *elem_map = new int[num_quads];
  int error = ex_get_map(exo_in, elem_map);
  if (error == 0) {
    printf("\n");
    int cnt = 0;
    int j;
    for (j = 0; j < num_quads; j++) {
      if (!verbose && j > 99)
        break;
          
      if (cnt == 0) {
        printf("\n%5d:    ", j + 1);
        cnt = 10;
      }
      printf(" %5d", elem_map[j]);
      --cnt;
    }
    printf("\n");
  }
  
  if (error == 0 && update)
    error = ex_put_map(exo_out, elem_map);
  
    // clean up
  delete [] elem_map;
  
  return error;
}

void print_property(int num_props, char* prop_name, int* prop_values)
{
  printf("\t%s:", prop_name);
  int i;
  for (i = 0; i < num_props; i++)
    printf(" %d", prop_values[i]);
  printf("\n");
}
