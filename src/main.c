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
    Main function.
*/
int main(int argc, char **argv)
{
    //this populates char output_dir[]
    if(!parse_options(argc, argv))
        return EXIT_FAILURE;
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
        // update_tangle contains calculation of normals, binormals and tangents.
        // Therefore we need to enforce boundaries, remeh (to keep right distance between points) 
        // and eliminate small loops (5 points needed for updating normals).
        reconnect(tangle, 0.0, rec_dist, reconnection_angle_cutoff); // this also enfrces boundaries, put it into separate function
        eliminate_small_loops(tangle, small_loop_cutoff);
        enforce_boundaries(tangle);
        remesh(tangle, global_dl_min, global_dl_max);
        update_tangle(tangle, 0);
        save_tangle(filename, tangle);
		fflush(stdout);
		sprintf(filename, "%s/frame0000.dat", output_dir);
		save_tangle(filename, tangle);
		fflush(stdout);
    }

    // define variables about tangle
    int shot = frame_shot - 1;
    int frame = 1;
    int recs = 0;
    int nrec = 0;
    int Np = number_of_used_nodes(tangle);

    // define time variables
    struct timespec t0, t_new, t_old;
    clockid_t clock = CLOCK_MONOTONIC;
    clock_gettime(clock, &t0);
    clock_gettime(clock, &t_new);
    clock_gettime(clock, &t_old);
    double time = 0;
    double delta_t1, delta_t2;

    // main siulation loop
    for(int k=0; Np > 0 && frame < frame_shot_total; ++k)
    {
        inject_vortices(tangle, time);
        nrec = reconnect(tangle, time, rec_dist, reconnection_angle_cutoff);
        recs += nrec;
        eliminate_small_loops(tangle, small_loop_cutoff);
        smooth(tangle);
        enforce_boundaries(tangle);
        remesh(tangle, global_dl_min, global_dl_max); // smooth has to be beforte update_tangle to remove points near each other
        update_tangle(tangle, time); // This needs to have all points inside the box because of velocities due to boundary images.
        rk4_step(tangle, time, global_dt);
        enforce_boundaries(tangle); // before saving and before recconection

        // Save frame.
        if(!shot)
	    {
            // Save tangle.
            sprintf(filename, "%s/frame%04d.dat", output_dir, frame);
            save_tangle(filename, tangle);
            frame++;
            shot = frame_shot;

            // Estimation of time left.
            clock_gettime(clock, &t_new);
            delta_t1 = time_diff(&t_old, &t_new);
            delta_t2 = time_diff(&t0, &t_new);
            clock_gettime(clock, &t_old);
            printf("Step %d, time = %g, recs: %d, Np: %d\n", k, time, recs, Np);
            printf("Current time per shot = %f s, avarage time per shot = %f s,avarage time left = %f min\n", delta_t1, delta_t2/frame, delta_t1*(frame_shot_total-frame)/60);
        }

        // Update number of points, number of steps till next shot, time.
        Np = number_of_used_nodes(tangle);
        shot--;
        time += global_dt;

        // Flush stdout.
        fflush(stdout);
    }

    // Write elapsed time per simulation.
    clock_gettime(clock, &t_new);
    printf("Elapsed seconds: %f\n", time_diff(&t0, &t_new));
    
    // Free the memory.
    free_tangle(tangle);
    free(tangle);

    return EXIT_SUCCESS;
}
