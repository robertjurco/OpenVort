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

void add_random_loops_centered(struct tangle_state *tangle, int N, double avg_radius, double variation, double dist);
void add_random_loops(struct tangle_state *tangle, int N, double avg_radius, double variation);

// line loops

void add_line_in_zdir(struct tangle_state *tangle, double x, double y, int direction);
void add_line(struct tangle_state *tangle, const struct segment *seg, int pin_wall_1, int pin_wall_2);
void add_line_with_oscilation(struct tangle_state *tangle, const struct segment *seg, int pin_wall_1, int pin_wall_2, double amplitude, double wavelength);
void add_line_with_two_oscilations(struct tangle_state *tangle, const struct segment *seg, int pin_wall_1, int pin_wall_2, double amplitude_1, double wavelength_1, double amplitude_2, double wavelength_2);
void add_line_in_zdir_with_two_oscilations_with_periodic_boundary(struct tangle_state *tangle, double x, double y, int direction, double amplitude_1, int number_of_waves_1, double amplitude_2, int number_of_waves_2);

// two random straight lines with opposite direction

void random_straight_lines_in_zdir(struct tangle_state *tangle, int N);

// inject vorticies

void inject_vortices(struct tangle_state *tangle, double t);
void inject_loop(struct tangle_state *tangle, double t, double frequency);
void inject_line_pairs(struct tangle_state *tangle, double t, double frequency, int line_injection_n, int polarized);

#endif //VORTEX_INJECTION_H
