#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <math.h>

#include "tangle.h"

extern int restart;
extern char restart_path[];

extern char output_dir[];

int parse_options(int argc, char **argv);
int setup_outdir(const char *dirname);

void save_tangle(const char *filename, struct tangle_state *tangle);
int load_tangle(const char *filename, struct tangle_state *tangle);

void print_usage(const char *prog_name);

#endif //UTILS_H
