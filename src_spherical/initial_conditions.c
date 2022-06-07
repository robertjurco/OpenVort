#include "tangle.h"
#include "vortex_constants.h"
#include "vortex_injection.h"
#include "vortex_smooth.h"

#include <stdio.h>

//declared as extern in vortex_constants.h

// Quantum properties.
double VORTEX_WIDTH = 1e-8; // Width of the vortex [cm].
double KAPPA = 9.97e-4; // Kappa [cm^2/s].

// General properties.
int frame_shot = 50; // How often to save a snapshot of the tangle.
int frame_shot_total = 3000; // Total number of saved snapshots (frames).
int global_num_threads = 8; // Number of threads of processor to use.

// Load tangle from init.dat located in output directory.
int load_tangle_from_init = 0;
int number_of_zeroth_frame = 0;

// For T=1.5K : alpha=0.72 ; alpha_p=0.01766 ; rho_n=16.17 ; rho_s = 129;
// For T=1.3K : alpha=0.034 ; alpha_p=1.383e-2 ; rho_n=6.522 ; rho_s = 138.60;

// Mutual friction.
int use_mutual_friction = 1; // Use mutual friction.
double alpha = 0.034; // Alpha coefficient of mutual friction.
double alpha_p = 1.383e-2; // Alpha prime coefficient of mutual friction.

// Densities [g/cm^3].
double rho_n = 6.522; // Normal fluid.
double rho_s = 138.60; // Superfluid.

// Resolution of the simulation.
double global_dt = 1e-4; // Time per one step, usually between 2e-5 s and 1e-4 s.
double global_dl_min = 0.5e-3; // Minimal length of discretisation, usually 0.5e-3.
double global_dl_max = 1.5e-3; // Maximal length of discretisation, usually 1.5e-3. This is at least twice the minimal distance.

// Reconnections and small loop cutoff.
double reconnection_angle_cutoff = DEG2RAD(5); // Minimal angle for reconnection, in radians, usually 5 degrees.
double rec_dist = 4e-4; // Recconection distance. Everythng closer then rec_dist is recconected. Has to be less then dl_min but more then dl_min/sqrt(2).
double rec_dist_with_walls = 2e-4; // Recconection distance with walls.
int small_loop_cutoff = 5; // Minimal number of points in the loop. If it reaches this number (usualy 5), the loop is deleted.

// Elimination of loops near the origin (for spherical flows).
int eliminate_origin_loops = 0; // Eliminate loops near origin (default off). The cutoff distance under which to eliminate loops is eliminate_loops_origin_cutoff.
double eliminate_loops_origin_cutoff = 3e-2; // The cutoff distance from origin under which to eliminate loops in spherical flow. About 3e-2 cm.

// Elimination of loops near the z axis (for cylindrical flows).
int eliminate_zaxis_loops = 0; // Eliminate loops near zaxis (default off). The cutoff distance under which to eliminate loops is eliminate_loops_zaxis_cutoff.
double eliminate_loops_zaxis_cutoff = 2e-2; // The cutoff distance from zaxis under which to eliminate loops in cylindrial flow. About 2e-2 cm.

// Pinning mode (PINNED or PINNED_SLIP).
int pin_mode = PINNED;

// Frequency filter.
int use_freq_filter = 1; // Bool, use frequency filter or not.
double freq_to_cutoff = 4.0/5; // Percentage of frequencies to cutoff (do not forget that we need decimal for double result).

// Barnes-Hut tree approximation.
double BH_resolution = 0.4; // Resolution of Barnes-Hut tree approximation, default 0.4.

// Power output of the heat source.
double q_dot = 0.03; // Power output of the heat source. // 0.05
double temp_bath = 1.65; // Temperature at heater.
int use_temperature = 1; // Use heat source and temperature distribution to caculate densities and mutual frictions.
int full_sphere = 1; // Are we simulating full sphere?

// Outside boundary.
int outer_surface = 0; // Outer surface is wall (1), open (0). If open everything fully behind it is killed.

//whether and from what file are we restarting the calculation
#define PATH_LEN 256
int restart = 0;
char restart_path[PATH_LEN] = "";

void setup_init_conditions(struct tangle_state *tangle)
{
	// Checks if the rec_dist has right value.
	if (rec_dist > global_dl_min || rec_dist < global_dl_min/sqrt(2))
	{
		printf("Wrong recconection distance. Recconection distance has to be more then dl_min/sqrt(2) but less then dl_min.");
		exit(EXIT_FAILURE);
	}
	// Checks if the global_dl has right value.
	if (global_dl_max < 2 * global_dl_min)
	{
		printf("Wrong discretisation distance. dl_max has to be at least twice as big as dl_min.");
		exit(EXIT_FAILURE);
	}
	
    //setup domain size
	tangle->domain_section.inner_radius = 0.1;
	tangle->domain_section.outer_radius = 0.5;
	tangle->domain_section.azimut_angle = DEG2RAD(90);
	tangle->domain_section.polar_angle = DEG2RAD(180);
	//full spehere? Put there 90, 180 and turn of image tangles

    //add vortices
    if (!load_tangle_from_init) 
	{
		add_random_loops(tangle, 40, 0.04, 0.7);
    }
}
