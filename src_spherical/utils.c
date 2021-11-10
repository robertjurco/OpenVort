#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>

#include "utils.h"
#include "tangle.h"

#define PATH_LEN 256

//+1 for the \0
char output_dir[PATH_LEN+1];

struct option options[] = {
    {"output", required_argument, NULL, 'o'},
    {0, 0, 0, 0}
};

int parse_options(int argc, char **argv)
{
    int c;
    int have_outdir = 0;

    c = getopt_long(argc, argv, "o:", options, NULL);
    
    if (c == 'o') {
        have_outdir = 1;
        strncpy(output_dir, optarg, PATH_LEN);
    }
    else goto failure;


    if(!have_outdir)
        goto failure;

    //trim the trailing / off of output_dir, if it's there
    char *x = output_dir;
    while(*(x+1)) x++;
    if(*x == '/') *x = 0;

    return 1;

    failure:
        print_usage(argv[0]);
        return 0;
}

int setup_outdir(const char *dirname)
{
    DIR *test = opendir(dirname);

    if(test){
        //the directory exist and we can open it
        closedir(test);
        return 1;
    } else {
        //either it doesn't exist or we can't open it
        //try creating it
        if(mkdir(dirname, 0777) == -1)
            return 0;
    }

    return 1;
}

void write_vector(FILE *stream, struct vec3 *v)
{
  fprintf(stream, "%.15g\t%.15g\t%.15g",
	  v->p[0], v->p[1], v->p[2]);
}

void save_point(FILE *stream, int vort_idx, const struct tangle_state *tangle, int i)
{
    fprintf(stream, "%d\t", vort_idx);
    write_vector(stream, tangle->vnodes + i);
    fprintf(stream, "\t");
    write_vector(stream, tangle->vels + i);
    fprintf(stream, "\t");
    write_vector(stream, tangle->tangents + i);
    fprintf(stream, "\t");
    write_vector(stream, tangle->normals + i);
    fprintf(stream, "\t%d\t%d\t%d\t%d\t%d", i,
        tangle->connections[i].reverse,
        tangle->connections[i].forward,
        tangle->status[i].status,
        tangle->status[i].pin_surface);
    fprintf(stream, "\n");
}

int load_point(FILE *stream, struct tangle_state *tangle)
{
    int zero;
    double pos0, pos1, pos2;
    double vel0, vel1, vel2;
    double tangent0, tangent1, tangent2;
    double normal0, normal1, normal2;
    int id;
    int reverse;
    int forward;
    int status;
    int pin_surface;

    int ret = fscanf(stream, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%d\t%d\t%d\t%d\t%d",
        &zero,
        &pos0,      &pos1,     &pos2,
 	    &vel0,      &vel1,     &vel2,
 	    &tangent0,  &tangent1, &tangent2,
 	    &normal0,   &normal1,  &normal2,
        &id,        &reverse,  &forward,
        &status,    &pin_surface);

    if(ret == EOF) return -1;

    while (id > tangle->N - 1) {
        expand_tangle(tangle, 2*tangle->N);
    }

    tangle->vnodes[id] = vec3(pos0, pos1, pos2);
    tangle->vels[id] = vec3(vel0, vel1, vel2);
    tangle->tangents[id] = vec3(tangent0, tangent1, tangent2);
    tangle->normals[id] = vec3(normal0, normal1, normal2);
    tangle->connections[id].reverse = reverse;
    tangle->connections[id].forward = forward;
    tangle->status[id].status = status;
    tangle->status[id].pin_surface = pin_surface;

  return id;
}

void save_tangle(const char *filename, struct tangle_state *tangle)
{
  int vortex_idx = 0;
  int *visited = (int *)calloc(tangle->N, sizeof(int));
  FILE *stream = fopen(filename, "w");

  if(!stream) {
      fprintf(stderr, "Can't open file %s\n", filename);
      return;
  }

  for(int k=0; k < tangle->N; ++k)
  {
       if(!visited[k])
	     {
	          if(tangle->status[k].status == EMPTY)
	          {
	              visited[k] = 1;
	              continue;
	          }

	          int first = k;
	          int curr  = k;
	          while(tangle->connections[curr].forward != first)
	          {
	               save_point(stream, vortex_idx, tangle, curr);
	               visited[curr] = 1;
	               curr = tangle->connections[curr].forward;
	               if(curr < 0)
		             {
		                   curr = tangle->connections[first].reverse;
		                   while(curr > 0 && tangle->connections[curr].reverse > 0)
		                   {
		                         save_point(stream, vortex_idx, tangle, curr);
		                         visited[curr] = 1;
		                         curr = tangle->connections[curr].reverse;
		                   }
		                   break;
		             }
	          }
	          if(curr > 0)
	          {
	               save_point(stream, vortex_idx, tangle, curr);
	               visited[curr] = 1;
	          }
	          vortex_idx++;
	     }
  }

  free(visited);
  fclose(stream);
}

int load_tangle(const char *filename, struct tangle_state *tangle)
{
    for (int i = 0; i < tangle->N; i++) {
        tangle->connections[i].reverse = -1;
        tangle->connections[i].forward = -1;
        tangle->status[i].status = EMPTY;
        tangle->status[i].pin_surface = NOT_A_SURFACE;
    }

    FILE *file = fopen(filename, "r");
    while(load_point(file, tangle) >=0 );
    fclose(file);

    return 0;
}

void print_usage(const char *prog_name)
{
    printf("Usage:\n");
    fputs(prog_name, stdout);
    printf(" -o | --output <output directory>\n");
}
