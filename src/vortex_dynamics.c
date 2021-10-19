#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "vortex_dynamics.h"
#include "vec3_math.h"
#include "tangle.h"
#include "vortex_constants.h"
#include "octree.h"

#include "boundaries.h"

/*
	Update the simulation using Euler 2.
	@param result: The updated tangle.
	@param tangle: The original tangle.
	@param dt: The time step.
	@param use_vels: Velocity array for nodes to use instead on velocities saved in tangle. This is used in Runge Kutta 4 method.
*/
void euler_step2(struct tangle_state *result, const struct tangle_state *tangle, double dt, const struct vec3 *use_vel)
{
	int k;
	
	for(k=0; k<tangle->N; ++k)
	{
		// Skip empty nodes.
		if(tangle->status[k].status == EMPTY) continue;

		// Update.
		struct vec3 move;
		if(!use_vel)
		{
	  		vec3_mul(&move, &tangle->vels[k], dt);
	  		vec3_add(&result->vnodes[k], &tangle->vnodes[k], &move);
		}
		else
		{
	  		vec3_mul(&move, &use_vel[k], dt);
	  		vec3_add(&result->vnodes[k], &tangle->vnodes[k], &move);
		}
	}
}

/*
	Update the simulation using Runge-Kutta 4.
	@param result: The tangle to update.
	@param t: Current time of the simulation.
	@param dt: The time step.
*/
void rk4_step(struct tangle_state *result, double t, double dt)
{
	// Make constant copy of the tangle.
    const struct tangle_state *tangle;
    tangle = result;

	// Number of nodes.
    int N;
	N = tangle->N;

	// Create tangles for Runge-Kutta 4: k2 through k4, we already have k1 in tangle.
    struct tangle_state rk_state[3];
    for(int k=0; k < 3; ++k)
    {
        create_tangle(&rk_state[k], N);
        // Only copy the stuff we actualy need.
        rk_state[k].bimg = tangle->bimg;
        rk_state[k].box = tangle->box;
        for(int kk = 0; kk < N; ++kk)
  		{
  	  		rk_state[k].status[kk] = tangle->status[kk];
  	  		rk_state[k].connections[kk] = tangle->connections[kk];
  		}
    }

    // Calculate k2 from k1.
    euler_step2(&rk_state[0], tangle, dt/2, NULL);
    update_tangents_normals(&rk_state[0]);
    update_velocities(&rk_state[0], t + dt/2);

    if (rk_state[0].N > N) printf("ERROR: rk_state[0]->N is greater then tangle->N\n");

    // Calculate k3 from k2..
    euler_step2(&rk_state[1], tangle, dt/2, rk_state[0].vels);
    update_tangents_normals(&rk_state[1]);
    update_velocities(&rk_state[1], t+dt/2);

    if (rk_state[1].N > N) printf("ERROR: rk_state[1]->N is greater then tangle->N\n");

    // Calculate k4 from k3.
    euler_step2(&rk_state[2], tangle, dt/2, rk_state[1].vels);
    update_tangents_normals(&rk_state[2]);
    update_velocities(&rk_state[2], t + dt);

	// do RK4
    for(int k = 0; k < N; ++k)
    {
		// Skip empty nodes.
        if(tangle->status[k].status == EMPTY) continue;

        // Move contains the full step, a are partial steps.
        struct vec3 move, a;
        vec3_mul(&move, &tangle->vels[k], dt/6);  // k1

        vec3_mul(&a, &rk_state[0].vels[k], dt/3); // k2
        vec3_add(&move, &move, &a);

        vec3_mul(&a, &rk_state[1].vels[k], dt/3); // k3
        vec3_add(&move, &move, &a);

        vec3_mul(&a, &rk_state[2].vels[k], dt/6); // k4
        vec3_add(&move, &move, &a);

		// Save the result.
        vec3_add(&result->vnodes[k], &result->vnodes[k], &move);
    }

	// Free the rk_state tangles.
    for(int k = 0; k < 3; ++k) free_tangle(&rk_state[k]);
}

/*
	Reconnect vortices at k-th and l-th nodes with each other.
	@param tangle : Tangle structure that holds all information about the nodes.
	@param k: k-th node.
	@param l: l-th node.
	@returns Return 1 if sucesfully reconnected, 0 otherwise.
*/
int do_reconnection(struct tangle_state* tangle, size_t k, size_t l)
{
	size_t kf, kr, lf, lr;

	double cf1, cf2;

	kf = tangle->connections[k].forward;
	kr = tangle->connections[k].reverse;

	lf = tangle->connections[l].forward;
	lr = tangle->connections[l].reverse;

	#define PSEG(x,y) seg_pwrap(tangle->vnodes+x, tangle->vnodes+y, &tangle->box)
	struct segment seg_kl = PSEG(k, l);
	struct segment seg_kflr = PSEG(kf, lr);
	struct segment seg_kkf = PSEG(k, kf);
	struct segment seg_llr = PSEG(l, lr);

	struct segment seg_krlf = PSEG(kr, lf);
	struct segment seg_kkr = PSEG(k, kr);
	struct segment seg_llf = PSEG(l, lf);
	#undef PSEG

	cf1 = segment_len(&seg_kl) + segment_len(&seg_kflr) - segment_len(&seg_kkf) - segment_len(&seg_llr);

	cf2 = segment_len(&seg_kl) + segment_len(&seg_krlf) - segment_len(&seg_kkr) - segment_len(&seg_llf);

	if (cf1 > 0 && cf2 > 0) return 0; //the reconnection increases the size, don't do it;

	if (cf1 < cf2)
	{
		tangle->connections[k].forward = l;
		tangle->connections[l].reverse = k;

		tangle->connections[kf].reverse = lr;
		tangle->connections[lr].forward = kf;
	}
	else
	{
		tangle->connections[l].forward = k;
		tangle->connections[k].reverse = l;

		tangle->connections[lf].reverse = kr;
		tangle->connections[kr].forward = lf;
	}

	return 1;
}

/*
	Run through all the pairs of nodes and check their distance if they are not
	immediate neighbors. Reconnect them if they are close enough. 
	This could, and should, in the future include some smarter estimation of possibility
	of a reconnections -- i.e., where in the BSP tree the two nodes are and only check them
	if can, in principle, be close enough.
	@param tangle: Tangle structure holding all the information about the vorices.
	@param t: Current time of the simulation.
 */
int reconnect(struct tangle_state *tangle, double t)
{
    int k, l;
    struct vec3 *v1, *v2; //points under test
    struct vec3 d1, d2; //direction vectors from v1, v2
    double calpha; //cosine of alpha between d1, d2
    int Nrecs = 0;

    Nrecs = 0;
    for(k = 0; k < tangle->N; ++k) tangle->recalculate[k] = 0;

    for(k = 0; k<tangle->N; ++k)
    {
        //do not reconnect points attached to the walls

        if(tangle->status[k].status == EMPTY 		||
            tangle->status[k].status == PINNED 		||
	 		tangle->status[k].status == PINNED_SLIP ||
            tangle->recalculate[k])
				continue;

		//skip empty nodes and nodes that reconnected in this pass
    	for(l=k+1; l<tangle->N; ++l)
		{
            if(tangle->status[l].status == EMPTY        ||
                tangle->status[l].status == PINNED      ||
                tangle->status[l].status == PINNED_SLIP ||
                tangle->connections[k].forward == l     ||
                tangle->connections[k].reverse == l     ||
                tangle->recalculate[l])
					continue;

	  		struct segment seg = seg_pwrap(tangle->vnodes + k, tangle->vnodes + l, &tangle->box);
	  		v1 = &seg.r1;
	  		v2 = &seg.r2;

	  		if(vec3_dist(v1, v2) > rec_dist) continue;
	  		//the nodes are close and they are not neighbors
	  		//now check the angle

	  		update_tangent_normal(tangle, k);
	  		update_tangent_normal(tangle, l);

	  		d1 = tangle->tangents[k];
	  		d2 = tangle->tangents[l];

	  		//In some cases vortex-vortex reconnections can create tiny loops (n ~ 3)
	  		//pinned on the wall which have colinear points for which calculation of the tangents fails.
	  		//These result in tangents of 0 length. If the tangent (which should be 1) is small
	  		//we just ignore this point because this loop will be removed on the next sweep of eliminate_small_loops
	  		//The 0.5 threshold is arbitrary. It will always be roughly 1 or roughly 0.
	  		if(vec3_len(&d1) < 0.5 || vec3_len(&d2) < 0.5) continue;
	  		//normalized dot -- just the cosine of the angle between vectors
	  		calpha = vec3_ndot(&d1, &d2);

	  		//We want to reconnect for angles LARGER than rec_angle
	  		//(i.e., parallel lines do not reconnect) and cosine is
	  		//a decreasing function.
	  		if(calpha > cos(reconnection_angle_cutoff))
	    		continue; //angle too small

	  		//Do not reconnect if the nodes are getting further apart from each other
	  		struct vec3 dx, dv;
	  		//we only update velocities if we really need them
	  		update_velocity(tangle, k, t, NULL);
	  		update_velocity(tangle, l, t, NULL);
	  		vec3_sub(&dx, v1, v2);
	  		vec3_sub(&dv, &tangle->vels[k], &tangle->vels[l]);

	  		if(vec3_dot(&dx, &dv) > 0)
				continue; //the points are getting further from each other

	  		//save the old connections because we might need them later
	  		int knext = tangle->connections[k].forward;
	  		int kprev = tangle->connections[k].reverse;
	  		int lnext = tangle->connections[l].forward;
	  		int lprev = tangle->connections[l].reverse;

	  		//angle is not too close and we can finally reconnect points k and l
	  		if(!do_reconnection(tangle, k, l))
	    		continue; //do_reconnect additionally checks if the total length doesn't increase

	  		//flag the neighbourhood as tainted so that it doesn't flash
	  		//back and forth in a single pass
	  		tangle->recalculate[k]++;
	  		tangle->recalculate[tangle->connections[k].forward]++;
	  		tangle->recalculate[tangle->connections[k].reverse]++;
	  		tangle->recalculate[knext]++;
	  		tangle->recalculate[kprev]++;

	  		tangle->recalculate[l]++;
	  		tangle->recalculate[tangle->connections[l].forward]++;
	  		tangle->recalculate[tangle->connections[l].reverse]++;
	  		tangle->recalculate[lnext]++;
	  		tangle->recalculate[lprev]++;

	  		Nrecs++;

	  		break;
		}
    }
    return Nrecs;
}

