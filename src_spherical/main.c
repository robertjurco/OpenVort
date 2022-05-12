#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#include <fenv.h>

#include "vec3_math.h"
#include "tangle.h"
#include "vortex_dynamics.h"
#include "vortex_constants.h"
#include "initial_conditions.h"
#include "utils.h"
#include "vortex_injection.h"
#include "vortex_smooth.h"

#include "temperature_functions.h"

/**
    Returns the time difference between two speified times in nanoseconds.
    @param t0 The old time.
    @param t1 The new time.
    @return The time difference in nanoseconds.
*/
double time_diff(struct timespec *t0, struct timespec *t1)
{
    double ns0 = t0->tv_sec*1e9 + t0->tv_nsec;
    double ns1 = t1->tv_sec*1e9 + t1->tv_nsec;

    return (ns1 - ns0)/1e9;
}

/**
	Save one frame of the simulation.
	@param tangle: Tangle holding all informations about vortices we want to save.
	@param step: Current step of the simulation.
	@param string: String to add at the end of the file.
	@param t_0: The time whe simulation started.
	@param t_new: The time of the curret step.
	@param t_old: The time of last step.
	@return The time difference in nanoseconds.
*/
void make_shot(struct tangle_state *tangle, int step, char* string, struct timespec *t_0, struct timespec *t_new, struct timespec* t_old)
{
	int frame = step / frame_shot;

	// Is it right time to make a shot?
	if (step % frame_shot != 0) return;

	// Save tangle.
	char filename[128];
	sprintf(filename, "%s/frame%04d%s.dat", output_dir, frame + number_of_zeroth_frame, string);
	save_tangle(filename, tangle);

	// Estimation of time left.
	double delta_t1 = time_diff(t_old, t_new);
	double delta_t2 = time_diff(t_0, t_new);
	int Np = number_of_used_nodes(tangle);
	
	// In case of frame > 0.
	if (frame > 0)
	{
		printf("Step %d, time = %g, Np: %d\n", step, step * global_dt, Np);
		printf("Current time per shot = %f s, avarage time per shot = %f s,avarage time left = %f min\n", delta_t1, delta_t2 / frame, delta_t1 * (frame_shot_total - frame) / 60);
	}
}

void under_radii(struct tangle_state* tangle)
{
	for (int i = 0; i < tangle->N; i++)
	{
		if (tangle->status[i].status == EMPTY) continue;
		if (vec3_len(&tangle->vnodes[i]) < tangle->domain_section.inner_radius * 0.95) {
			printf("Hele je pod inner radius %g\n", vec3_len(&tangle->vnodes[i]));
		}
	}
}


/**
    Main function.
*/
int main(int argc, char **argv)
{
    // This populates char output_dir[].
    if(!parse_options(argc, argv)) return EXIT_FAILURE;
    setup_outdir(output_dir);
    char filename[128];

    srand48(time(NULL));
    feenableexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO);

    //setup tangle
    struct tangle_state *tangle = (struct tangle_state*)malloc(sizeof(struct tangle_state));
    create_tangle(tangle, 512); //512
    setup_init_conditions(tangle);
	
    //load tangle or save generated one
    fflush(stdout);
    sprintf(filename, "%s/init.dat", output_dir);
    if (load_tangle_from_init) {
        load_tangle(filename, tangle);
    } else {
		enforce_wall_boundaries(tangle); // Goes first in case of vortex initialized outside of the domain. 
		eliminate_small_loops(tangle); // Enforce wall boundaries may generate small loops on walls, this may interfere with recconect and broke tangents and normals.
		reconnect(tangle, 0.0);
		eliminate_small_loops(tangle);  // Reconect may generate small loops on walls, this may interfere with smooth and broke tangents and normals.

		smooth(tangle); // No small loops, everithng reconnected, time to smooth.
		enforce_periodic_boundaries(tangle); // Smooth destroyes periodicity.
		remesh(tangle, global_dl_min, global_dl_max); // Remesh has to go beforte update_tangle to remove points near each other.
		eliminate_small_loops(tangle);  // Remesh may generate small loops on walls, this may interfere with smooth and broke tangents and normals.

		update_tangle(tangle, 0.0); // This needs to have all points inside the box because of velocities due to boundary images.

		save_tangle(filename, tangle);
		fflush(stdout);
    }

    // define time variables
    struct timespec t_0, t_new, t_old;
    clockid_t clock = CLOCK_MONOTONIC;
    clock_gettime(clock, &t_0);
    clock_gettime(clock, &t_new);
    clock_gettime(clock, &t_old);
    double time = 0;

    // main siulation loop
    for(int k = 0; number_of_used_nodes(tangle) > 0 && k < frame_shot * frame_shot_total; ++k)
    {
		// Make a shot needs to go first as we are here starting with in ccase of loaded tangle.
		// init.dat and frame0000.dat are always the same.
		make_shot(tangle, k, "", &t_0, &t_new, &t_old);

        //inject_vortices(tangle, time);
		enforce_wall_boundaries(tangle); // Goes first in case of vortex initialized outside of the domain, also returns PINNED_SLIP points onto the wall.
		eliminate_small_loops(tangle); // Enforce wall boundaries may generate small loops on walls, this may interfere with recconect and broke tangents and normals.
		reconnect(tangle, 0.0);
		eliminate_small_loops(tangle);  // Reconect may generate small loops on walls, this may interfere with smooth and broke tangents and normals.

		smooth(tangle); // No small loops, everithng reconnected, time to smooth.
		enforce_periodic_boundaries(tangle); // Smooth destroyes periodicity.
		remesh(tangle, global_dl_min, global_dl_max); // Remesh has to go beforte update_tangle to remove points near each other.
		eliminate_small_loops(tangle);  // Remesh may generate small loops on walls, this may interfere with smooth and broke tangents and normals.

		update_tangle(tangle, time); // This needs to have all points inside the box because of velocities due to boundary images.

		rk4_step(tangle, time, global_dt);

		enforce_periodic_boundaries(tangle); // Before saving.

		// Update time.
		time += global_dt;
		t_old.tv_sec = t_new.tv_sec;
		t_old.tv_nsec = t_new.tv_nsec;
		clock_gettime(clock, &t_new);

		// Flush stdout.
		fflush(stdout);
    }

    // Write elapsed time per simulation.
    clock_gettime(clock, &t_new);
    printf("Elapsed seconds: %f\n", time_diff(&t_0, &t_new));
    
    // Free the memory.
    free_tangle(tangle);
    free(tangle);

    return EXIT_SUCCESS;
}
