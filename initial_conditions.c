#include "tangle.h"
#include "vortex_constants.h"
#include "external_velocity.h"
#include "vortex_injection.h"
#include "vortex_smooth.h"

#include <stdio.h>

//declared as extern in vortex_constants.h

// Quantum properties.
double VORTEX_WIDTH = 1e-8; // Width of the vortex [cm].
double KAPPA = 9.97e-4; // Kappa [cm^2/s].

// General properties.
int frame_shot = 50; // How often to save a snapshot of the tangle.
int frame_shot_total = 100000; // Total number of saved snapshots (frames).
int global_num_threads = 6; // Number of threads of processor to use.

// Load tangle from init.dat located in output directory.
int load_tangle_from_init = 0;

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
double global_dl_max = 1.5e-3; // Maximal length of discretisation, usually 1.5e-3.

// Reconnections and small loop cutoff.
double reconnection_angle_cutoff = DEG2RAD(5); // Minimal angle for reconnection, in radians, usually 5 degrees.
double rec_dist = 4e-4; // Recconection distance. Everythng closer then rec_dist is recconected. Has to be less then dl_min but more then dl_min/sqrt(2).
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
double freq_to_cutoff = 4/5; // Percentage of frequencies to cutoff.

// Injecting loops at upper z-plane.
int loop_injection = 0; // Bool, inject or not.
double loop_injection_frequency = 1; // Number of loops injected per second (frequency).

// Line pair injection.
int line_injection = 0; // Bool, inject or not.
int line_injection_n = 1; // How many pairs to inject.
double line_injection_frequency = 1; // Number of pairs injected per second (frequency).
int line_injection_polarized = 0; // Bool, whether the injection is polarized. Off by default.

// Barnes-Hut tree approximation.
int use_BH = 1; // 1 for use of Barnes-Hut tree approximation.
double BH_resolution = 0.4; // Resolution of Barnes-Hut tree approximation, default 0.4.

// Velocity of normal fluid.
struct ext_vel_param vn_conf = {
                                .type = SIMPLE,
                                .strength = 4.775292513, //v_ns/(1-rho_n/rho_s),
                                .direction = {{1.0,0,0}}
                                };

// Velocity of boundary.
struct ext_vel_param vb_conf = {
                                .type = NO_FLOW
                                };

// Velocity of superfluid.
struct ext_vel_param vs_conf = {
                                .type = SIMPLE,
                                .strength = 0.2247074875, //-v_ns/(1-rho_s/rho_n),
                                .direction = {{-1.0,0,0}}
                                };

// No external or boundary velocity.								
/*
struct ext_vel_param vn_conf = {
                                .type = NO_FLOW,
                                .strength = 0,
                                .cut_off = 0,
                                .direction = {{0,0,0}}
                                };
struct ext_vel_param vb_conf = {
                                .type = NO_FLOW,
                                .strength = 0,
                                .cut_off = 0,
                                .direction = {{0,0,0}}
                                };
struct ext_vel_param vs_conf = {
                                .type = NO_FLOW,
                                .strength = 0,
                                .cut_off = 0,
                                .direction = {{0,0,0}}
                                };
*/

//whether and from what file are we restarting the calculation
#define PATH_LEN 256
int restart = 0;
char restart_path[PATH_LEN] = "";

void set_walls_full(struct tangle_state *tangle, wall_type wall)
{
  for(int k=0; k<6; ++k)
        tangle->box.wall[k] = wall;
}

void setup_init_conditions(struct tangle_state *tangle)
{
	// Checks if the rec_dist has right value.
	if (rec_dist > global_dl_min || rec_dist < global_dl_min/sqrt(2))
	{
		printf("Wrong recconection distance. Recconection distance has to be more then dl_min/sqrt(2) but less then dl_min.");
		exit(EXIT_FAILURE);
	}
	
    //setup domain size
    tangle->box.bottom_left_back = vec3(-0.05,-0.05,-0.05);
    tangle->box.top_right_front = vec3(0.05,0.05,0.05);

    //setup boundariess

    /*
    tangle->bimg = periodic_6; //18, 26
    set_walls_full(tangle, WALL_PERIODIC);
    */

    
    tangle->bimg = wall_canal;
    set_walls_full(tangle, WALL_MIRROR);
    tangle->box.wall[RIGHT] = WALL_PERIODIC;
    tangle->box.wall[LEFT] = WALL_PERIODIC;
    

    /*
    tangle->bimg = wall_1_6;
    set_walls_full(tangle, WALL_PERIODIC);
    tangle->box.wall[DOWN] = WALL_MIRROR;
    */

    /*
    tangle->bimg = open_boundaries;
    set_walls_full(tangle, WALL_OPEN);
    */
	/*
    tangle->bimg = wall_2_2;
    set_walls_full(tangle, WALL_MIRROR);
    tangle->box.wall[LEFT] = WALL_PERIODIC;
    tangle->box.wall[RIGHT] = WALL_PERIODIC;
	*/
    //add vortices
    if (!load_tangle_from_init) {

    /*
    struct segment dir_seg1, dir_seg2;
    dir_seg1.r1 = vec3(0.0, 0.015, -0.1);
    dir_seg1.r2 = vec3(0.0, 0.015, 0.1);
    dir_seg2.r1 = vec3(0.0, -0.015, 0.1);
    dir_seg2.r2 = vec3(0.0, -0.015, -0.1);

    add_wall_dir_line(tangle, &dir_seg1);
    add_wall_dir_line(tangle, &dir_seg2);
    */

    /*
    add_line(tangle, 0, -0.0025, +1, 125);
    add_line(tangle, 0, 0.0025, +1, 125);
    */
    /*
    struct vec3 center1 = vec3(0,0,0);
    struct vec3 dir1 = vec3(-1,0,0);
    add_circle(tangle, &center1, &dir1, 0.05);
    */

    /*
    struct vec3 dir = vec3(-1,0,0);
    add_circle_oscilation(tangle, &dir, 0.05, 0.005, 0.02);
	*/

    /*

    struct segment dir_seg1;
    dir_seg1.r1 = vec3(0, 0, -0.1);
    dir_seg1.r2 = vec3(0, 0, 0.1);
    add_wall_dir_line_oscilation_xy(tangle, &dir_seg1, 0.02, 0.04); // k = 5
    */
 
   /*
    struct segment dir_seg1, dir_seg2, dir_seg3; // they are 0.06 appart from each other
    dir_seg1.r1 = vec3(0, 0.03464101615, -0.1);
    dir_seg1.r2 = vec3(0, 0.03464101615, 0.1);
    add_wall_dir_line_two_oscilations_xy(tangle, &dir_seg1, 0.004, 0.1, 0.004, 0.025); // k = 2 and 8
    dir_seg2.r1 = vec3(0.03, -0.01732050808, -0.1);
    dir_seg2.r2 = vec3(0.03, -0.01732050808, 0.1);
    add_wall_dir_line_two_oscilations_xy(tangle, &dir_seg2, 0.004, 0.06666666666, 0.004, 0.02857142857); // k = 3 and 7
    dir_seg3.r1 = vec3(-0.03, -0.01732050808, -0.1);
    dir_seg3.r2 = vec3(-0.03, -0.01732050808, 0.1);
    add_wall_dir_line_two_oscilations_xy(tangle, &dir_seg3, 0.004, 0.05, 0.004, 0.03333333333); // k = 4 and 6
    */

	struct vec3 center1 = vec3(0, 0.01, 0.05);
	struct vec3 dir1 = vec3(1, 0, 0);
	add_loop(tangle, &center1, &dir1, 0.04);

	struct vec3 center2 = vec3(0.03, 0.0, 0.05);
	struct vec3 dir2 = vec3(1, 0, 0);
	add_loop(tangle, &center2, &dir2, 0.02);

	struct vec3 center3 = vec3(0.02, 0.03, -0.05);
	struct vec3 dir3 = vec3(1, 0, 0);
	add_loop(tangle, &center3, &dir3, 0.02);

	struct vec3 center4 = vec3(0.04, -0.03, -0.05);
	struct vec3 dir4 = vec3(1, 0, 0);
	add_loop(tangle, &center4, &dir4, 0.015);

	
    /*
    add_random_loops(tangle, 20, 0.03, 0.25);
    */
    }
}
