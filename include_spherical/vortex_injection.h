/**
    Header for initialization and injecion of vorticies.
    @file vortex_injection.h
*/

#ifndef VORTEX_INJECTION_H
#define VORTEX_INJECTION_H

#include "tangle.h"

// one loop

void add_loop(struct tangle_state *tangle, struct vec3 *center, struct vec3 *dir, double r);
void add_loop_with_oscilation(struct tangle_state *tangle, struct vec3* center, struct vec3 *dir, double r, double amplitude, double wavelength);

// random loops

void add_random_loops(struct tangle_state* tangle, int N, double avg_radius, double variation);

#endif //VORTEX_INJECTION_H
