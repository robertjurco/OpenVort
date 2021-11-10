#ifndef VORTEX_CONSTANTS_H
#define VORTEX_CONSTANTS_H

#include <math.h>

#define DEG2RAD(X) ((X)*M_PI/180.0)

// These are defined with default values in vortex_configuration.c

// Quantum properties.
extern double VORTEX_WIDTH; // Width of the vortex [cm].
extern double KAPPA; // Kappa [cm^2/s].

// General properties.
extern int frame_shot; // How often to save a snapshot of the tangle.
extern int frame_shot_total; // Total number of saved snapshots (frames).
extern int global_num_threads; // Number of threads of processor to use.

// Load tangle from init.dat located in output directory.
extern int load_tangle_from_init;

// Mutual friction.
extern int use_mutual_friction; // Use mutual friction.
extern double alpha; // Alpha coefficient of mutual friction.
extern double alpha_p; // Alpha prime coefficient of mutual friction.

// Densities [g/cm^3].
extern double rho_n; // Normal fluid.
extern double rho_s; // Superfluid.

// Resolution of the simulation.
extern double global_dt; // Time per one step, usually between 2e-5 s and 1e-4 s.
extern double global_dl_min; // Minimal length of discretisation, usually 0.5e-3.
extern double global_dl_max; // Maximal length of discretisation, usually 1.5e-3.

// Reconnections and small loop cutoff.
extern double reconnection_angle_cutoff; // Minimal angle for reconnection, in radians, usually 5 degrees.
extern double rec_dist; // Recconection distance. Everythng closer then rec_dist is recconected. Has to be less then dl_min but more then dl_min/sqrt(2).
extern double rec_dist_with_walls; // Recconection distance with walls.
extern int small_loop_cutoff; // Minimal number of points in the loop. If it reaches this number (usualy 5), the loop is deleted.

// Pinning mode (PINNED or PINNED_SLIP).
extern int pin_mode;

// Frequency filter.
extern int use_freq_filter; // Bool, use frequency filter or not.
extern double freq_to_cutoff; // Percentage of frequencies to cutoff.

// Barnes-Hut tree approximation.
extern double BH_resolution; // Resolution of Barnes-Hut tree approximation, default 0.4.


#endif //VORTEX_CONSTANTS_H
